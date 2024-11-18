// MusicPlaylistManagerTest.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <limits>
#include <windows.h>
#include "LiveRecording.h"

using namespace std;

// Configuration constants
const int STATE_COUNT = 5;
const int MAX_SEQUENCE_LENGTH = 100;
const int OBSERVATION_SYMBOLS = 32;
const int WORD_COUNT = 5;
const int TRAINING_SET_SIZE = 30;
const int TEST_SET_SIZE = 10;
const int ITERATION_COUNT = 6;
const int BAUM_WELCH_ITERATIONS = 60;
const string DATA_DIRECTORY = "210101070_dataset/observation/";
const double MIN_PROBABILITY = 1e-15;

// Word mapping structure
struct WordMapping {
    const char* word;
    int index;
} wordMap[] = {
    {"create", 0},
    {"add", 1},
    {"play", 2},
    {"stop", 3},
	{"shuffle", 4}
};

// HMM model structure
struct HMMParameters
{
    double transitionMatrix[STATE_COUNT][STATE_COUNT];
    double emissionMatrix[STATE_COUNT][OBSERVATION_SYMBOLS];
    double initialProbs[STATE_COUNT];
};

struct ComputationMatrices
{
    int observationSeq[MAX_SEQUENCE_LENGTH];
    int seqLength;
    long double forwardProbs[MAX_SEQUENCE_LENGTH][STATE_COUNT];
    long double backwardProbs[MAX_SEQUENCE_LENGTH][STATE_COUNT];
    long double stateProbs[MAX_SEQUENCE_LENGTH][STATE_COUNT];
    long double transitionProbs[MAX_SEQUENCE_LENGTH][STATE_COUNT][STATE_COUNT];
};

ComputationMatrices compute;

// File handling functions
bool fileExists(const string &filepath)
{
    ifstream file(filepath);
    return file.good();
}

// Load matrices and vectors from files
void loadEmissionMatrixFromFile(const string &filepath, double mat[][OBSERVATION_SYMBOLS], int numRows, int numCols)
{
    ifstream inputFile(filepath.c_str());
    if (!inputFile)
    {
        cerr << "Error: Unable to open the file " << filepath.c_str() << endl;
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < numRows; ++i)
    {
        for (int j = 0; j < numCols; ++j)
        {
            inputFile >> mat[i][j];
        }
    }
    inputFile.close();
}

void loadTransitionMatrixFromFile(const string &filepath, double mat[][STATE_COUNT], int numRows, int numCols)
{
    ifstream inputFile(filepath.c_str());
    if (!inputFile)
    {
        cerr << "Error: Unable to open the file " << filepath.c_str() << endl;
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < numRows; ++i)
    {
        for (int j = 0; j < numCols; ++j)
        {
            inputFile >> mat[i][j];
        }
    }
    inputFile.close();
}

void loadVectorFromFile(const string &filepath, double vec[], int length)
{
    ifstream inputFile(filepath.c_str());
    if (!inputFile)
    {
        cerr << "Error: Unable to open the file " << filepath.c_str() << endl;
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < length; ++i)
    {
        inputFile >> vec[i];
    }
    inputFile.close();
}

void initializeModelFromFiles(HMMParameters &model)
{
    // Load matrices and vectors from files
    loadTransitionMatrixFromFile("files/A_ijMatrix.txt", model.transitionMatrix, 5, 5);
    loadEmissionMatrixFromFile("files/B_jkMatrix.txt", model.emissionMatrix, 5, 32);
    loadVectorFromFile("files/P_iVector.txt", model.initialProbs, 5);

    // Validate probabilities
    for (int i = 0; i < STATE_COUNT; i++)
    {
        double transitionSum = 0, emissionSum = 0;
        for (int j = 0; j < STATE_COUNT; j++)
        {
            transitionSum += model.transitionMatrix[i][j];
        }
        for (int j = 0; j < OBSERVATION_SYMBOLS; j++)
        {
            emissionSum += model.emissionMatrix[i][j];
        }
        if (abs(transitionSum - 1.0) > 1e-6 || abs(emissionSum - 1.0) > 1e-6)
        {
            cerr << "Warning: Row " << i << " probabilities do not sum to 1" << endl;
        }
    }
}

// Model averaging and evaluation functions
void averageModels(const vector<HMMParameters> &models, HMMParameters &result)
{
    int modelCount = (int)models.size();
    if (modelCount == 0)
        return;

    // Average transition matrix
    for (int i = 0; i < STATE_COUNT; i++)
    {
        for (int j = 0; j < STATE_COUNT; j++)
        {
            double sum = 0.0;
            for (int k = 0; k < modelCount; k++)
            {
                sum += models[k].transitionMatrix[i][j];
            }
            result.transitionMatrix[i][j] = sum / modelCount;
        }
    }

    // Average emission matrix
    for (int i = 0; i < STATE_COUNT; i++)
    {
        for (int j = 0; j < OBSERVATION_SYMBOLS; j++)
        {
            double sum = 0.0;
            for (int k = 0; k < modelCount; k++)
            {
                sum += models[k].emissionMatrix[i][j];
            }
            result.emissionMatrix[i][j] = sum / modelCount;
        }
    }

    // Average initial probabilities
    for (int i = 0; i < STATE_COUNT; i++)
    {
        double sum = 0.0;
        for (int k = 0; k < modelCount; k++)
        {
            sum += models[k].initialProbs[i];
        }
        result.initialProbs[i] = sum / modelCount;
    }
}

bool loadObservationSequence(const string &filepath)
{
    ifstream file(filepath);
    if (!file)
    {
        return false;
    }

    compute.seqLength = 0;
    int observation;
    while (file >> observation && compute.seqLength < MAX_SEQUENCE_LENGTH)
    {
        compute.observationSeq[compute.seqLength++] = observation;
    }
    file.close();
    return true;
}

// Core HMM functions
long double calculateSequenceProbability(const HMMParameters &model)
{
    // Initialize first time step
    for (int i = 0; i < STATE_COUNT; i++)
    {
        compute.forwardProbs[0][i] = model.initialProbs[i] *
                                     model.emissionMatrix[i][compute.observationSeq[0]];
    }

    // Forward recursion
    for (int t = 1; t < compute.seqLength; t++)
    {
        for (int j = 0; j < STATE_COUNT; j++)
        {
            compute.forwardProbs[t][j] = 0.0;
            for (int i = 0; i < STATE_COUNT; i++)
            {
                compute.forwardProbs[t][j] += compute.forwardProbs[t - 1][i] *
                                              model.transitionMatrix[i][j];
            }
            compute.forwardProbs[t][j] *= model.emissionMatrix[j][compute.observationSeq[t]];
        }
    }

    // Calculate total probability
    long double probability = 0.0;
    for (int i = 0; i < STATE_COUNT; i++)
    {
        probability += compute.forwardProbs[compute.seqLength - 1][i];
    }
    return probability;
}

// Evaluate sequence probability using a given model
long double evaluateSequence(const HMMParameters &model, const string &filepath)
{
    if (!loadObservationSequence(filepath))
    {
        return -numeric_limits<long double>::infinity();
    }
    return calculateSequenceProbability(model);
}

// Save and load model functions
void saveModelToFile(const HMMParameters &model, int digit)
{
    system("mkdir models 2> nul");

    ostringstream filepath;
    filepath << "models/command_" << wordMap[digit].word << ".txt";
    ofstream file(filepath.str());

    if (!file)
    {
        cerr << "Error: Could not save model for word " << wordMap[digit].word << "\n";
        return;
    }

    file << "Transition Matrix:\n";
    for (int i = 0; i < STATE_COUNT; i++)
    {
        for (int j = 0; j < STATE_COUNT; j++)
        {
            file << scientific << setprecision(10) << model.transitionMatrix[i][j] << " ";
        }
        file << "\n";
    }

    file << "\nEmission Matrix:\n";
    for (int i = 0; i < STATE_COUNT; i++)
    {
        for (int j = 0; j < OBSERVATION_SYMBOLS; j++)
        {
            file << scientific << setprecision(10) << model.emissionMatrix[i][j] << " ";
        }
        file << "\n";
    }

    file << "\nInitial Probabilities:\n";
    for (int i = 0; i < STATE_COUNT; i++)
    {
        file << scientific << setprecision(10) << model.initialProbs[i] << " ";
    }
    file << "\n";
    file.close();
}

bool loadModelFromFile(HMMParameters &model, int digit)
{
    ostringstream filepath;
    filepath << "models/command_" << wordMap[digit].word << ".txt";
    ifstream file(filepath.str());

    if (!file)
    {
        cerr << "Error: Could not load model for word " << wordMap[digit].word << "\n";
        return false;
    }

    string line;
    getline(file, line); // Skip header

    for (int i = 0; i < STATE_COUNT; i++)
    {
        for (int j = 0; j < STATE_COUNT; j++)
        {
            file >> model.transitionMatrix[i][j];
        }
    }

    getline(file, line); // Skip newline
    getline(file, line); // Skip header
    getline(file, line); // Skip newline

    for (int i = 0; i < STATE_COUNT; i++)
    {
        for (int j = 0; j < OBSERVATION_SYMBOLS; j++)
        {
            file >> model.emissionMatrix[i][j];
        }
    }

    getline(file, line); // Skip newline
    getline(file, line); // Skip header
    getline(file, line); // Skip newline

    for (int i = 0; i < STATE_COUNT; i++)
    {
        file >> model.initialProbs[i];
    }

    file.close();
    return true;
}

// Baum-Welch algorithm implementation
void baumWelch(HMMParameters &model, const string &filepath, int maxIterations = BAUM_WELCH_ITERATIONS)
{
    if (!loadObservationSequence(filepath))
    {
        return;
    }

    double previousProbability = 0.0;

    for (int iter = 0; iter < maxIterations; iter++)
    {
        // Forward algorithm
        long double currentProbability = calculateSequenceProbability(model);
        double probImprovement = fabs(currentProbability - previousProbability);

        if (iter > 0 && probImprovement == 0)
        {
            break;
        }
        previousProbability = currentProbability;

        // Backward algorithm
        for (int i = 0; i < STATE_COUNT; i++)
        {
            compute.backwardProbs[compute.seqLength - 1][i] = 1.0;
        }

        for (int t = compute.seqLength - 2; t >= 0; t--)
        {
            for (int i = 0; i < STATE_COUNT; i++)
            {
                compute.backwardProbs[t][i] = 0.0;
                for (int j = 0; j < STATE_COUNT; j++)
                {
                    compute.backwardProbs[t][i] += model.transitionMatrix[i][j] *
                                                   model.emissionMatrix[j][compute.observationSeq[t + 1]] *
                                                   compute.backwardProbs[t + 1][j];
                }
            }
        }

        // Compute Xi and Gamma
        for (int t = 0; t < compute.seqLength - 1; t++)
        {
            long double denominator = 0.0;
            for (int i = 0; i < STATE_COUNT; i++)
            {
                for (int j = 0; j < STATE_COUNT; j++)
                {
                    denominator += compute.forwardProbs[t][i] *
                                   model.transitionMatrix[i][j] *
                                   model.emissionMatrix[j][compute.observationSeq[t + 1]] *
                                   compute.backwardProbs[t + 1][j];
                }
            }

            for (int i = 0; i < STATE_COUNT; i++)
            {
                compute.stateProbs[t][i] = 0.0;
                for (int j = 0; j < STATE_COUNT; j++)
                {
                    if (denominator > 0)
                    {
                        compute.transitionProbs[t][i][j] = (compute.forwardProbs[t][i] *
                                                            model.transitionMatrix[i][j] *
                                                            model.emissionMatrix[j][compute.observationSeq[t + 1]] *
                                                            compute.backwardProbs[t + 1][j]) /
                                                           denominator;
                    }
                    else
                    {
                        compute.transitionProbs[t][i][j] = 0.0;
                    }
                    compute.stateProbs[t][i] += compute.transitionProbs[t][i][j];
                }
            }
        }

        // Handle final time step
        long double finalDenominator = 0.0;
        for (int i = 0; i < STATE_COUNT; i++)
        {
            finalDenominator += compute.forwardProbs[compute.seqLength - 1][i] *
                                compute.backwardProbs[compute.seqLength - 1][i];
        }

        for (int i = 0; i < STATE_COUNT; i++)
        {
            if (finalDenominator > 0)
            {
                compute.stateProbs[compute.seqLength - 1][i] =
                    (compute.forwardProbs[compute.seqLength - 1][i] *
                     compute.backwardProbs[compute.seqLength - 1][i]) /
                    finalDenominator;
            }
            else
            {
                compute.stateProbs[compute.seqLength - 1][i] = 0.0;
            }
        }

        // Re-estimate model parameters
        // Update Pi
        for (int i = 0; i < STATE_COUNT; i++)
        {
            model.initialProbs[i] = compute.stateProbs[0][i];
        }

        // Update A
        for (int i = 0; i < STATE_COUNT; i++)
        {
            for (int j = 0; j < STATE_COUNT; j++)
            {
                double numerator = 0.0, denominator = 0.0;
                for (int t = 0; t < compute.seqLength - 1; t++)
                {
                    numerator += compute.transitionProbs[t][i][j];
                    denominator += compute.stateProbs[t][i];
                }
                if (denominator > 0)
                {
                    model.transitionMatrix[i][j] = numerator / denominator;
                }
            }
        }

        // Update B
        for (int i = 0; i < STATE_COUNT; i++)
        {
            for (int k = 0; k < OBSERVATION_SYMBOLS; k++)
            {
                double numerator = 0.0, denominator = 0.0;
                for (int t = 0; t < compute.seqLength; t++)
                {
                    if (compute.observationSeq[t] == k)
                    {
                        numerator += compute.stateProbs[t][i];
                    }
                    denominator += compute.stateProbs[t][i];
                }
                if (denominator > 0)
                {
                    model.emissionMatrix[i][k] = numerator / denominator;
                    if (model.emissionMatrix[i][k] < MIN_PROBABILITY)
                    {
                        model.emissionMatrix[i][k] = MIN_PROBABILITY;
                    }
                }
            }

            // Normalize emission probabilities
            double rowSum = 0.0;
            for (int k = 0; k < OBSERVATION_SYMBOLS; k++)
            {
                rowSum += model.emissionMatrix[i][k];
            }
            if (rowSum > 1.0)
            {
                int maxIndex = 0;
                for (int k = 1; k < OBSERVATION_SYMBOLS; k++)
                {
                    if (model.emissionMatrix[i][k] > model.emissionMatrix[i][maxIndex])
                    {
                        maxIndex = k;
                    }
                }
                model.emissionMatrix[i][maxIndex] -= (rowSum - 1.0);
                if (model.emissionMatrix[i][maxIndex] < MIN_PROBABILITY)
                {
                    model.emissionMatrix[i][maxIndex] = MIN_PROBABILITY;
                }
            }
        }
    }
}

// Model training and testing functions
// Training function
void trainModels(vector<HMMParameters> &wordModels)
{
    cout << "Starting model training for all words...\n";

    for (int wordIdx = 0; wordIdx < WORD_COUNT; wordIdx++)
    {
        cout << "\nTraining word '" << wordMap[wordIdx].word << "'...\n";
        initializeModelFromFiles(wordModels[wordIdx]);

        for (int iteration = 0; iteration < ITERATION_COUNT; iteration++)
        {
            cout << "Global iteration " << (iteration + 1) << "/" << ITERATION_COUNT << "\n";
            vector<HMMParameters> iterationModels;

            for (int sample = 1; sample <= TRAINING_SET_SIZE; sample++)
            {
                ostringstream filepath;
                filepath << DATA_DIRECTORY << "210101070_W_" << wordIdx << "_" << sample << ".txt";

                HMMParameters currentModel = wordModels[wordIdx];
                baumWelch(currentModel, filepath.str());
                iterationModels.push_back(currentModel);
            }

            if (!iterationModels.empty())
            {
                averageModels(iterationModels, wordModels[wordIdx]);
            }
        }

        cout << "Saving model for word '" << wordMap[wordIdx].word << "'...\n";
        saveModelToFile(wordModels[wordIdx], wordIdx);
    }

    cout << "\nTraining and saving completed successfully.\n";
}

// Testing function
void testModels(const vector<HMMParameters> &wordModels)
{
    int correctCount = 0;
    int uncertainCount = 0;
    int totalTests = WORD_COUNT * TEST_SET_SIZE;

    cout << "\nStarting model testing phase...\n";

    for (int actualWord = 0; actualWord < WORD_COUNT; actualWord++)
    {
        int wordCorrect = 0;
        int wordUncertain = 0;

        cout << "\nTesting word '" << wordMap[actualWord].word << "'...\n";

        for (int sample = 1; sample <= TEST_SET_SIZE; sample++)
        {
            ostringstream filepath;
            filepath << DATA_DIRECTORY << "210101070_W_" << actualWord << "_"
                     << (TRAINING_SET_SIZE + sample) << ".txt";

            // Get probabilities for all models
            vector<pair<int, long double>> modelProbabilities;
            for (int modelWord = 0; modelWord < WORD_COUNT; modelWord++)
            {
                long double probability = evaluateSequence(wordModels[modelWord], filepath.str());
                modelProbabilities.push_back(make_pair(modelWord, probability));
            }

            // Find the best model
            auto bestModel = max_element(modelProbabilities.begin(), modelProbabilities.end(),
                                         [](const pair<int, long double> &a, const pair<int, long double> &b)
                                         { return a.second < b.second; });

            int predictedWord = bestModel->first;
            long double bestProb = bestModel->second;

            // Find second best for uncertainty check
            auto secondBest = max_element(modelProbabilities.begin(), modelProbabilities.end(),
                                          [bestModel](const pair<int, long double> &a, const pair<int, long double> &b)
                                          {
                                              return (a.first == bestModel->first || a.second < b.second) &&
                                                     (b.first != bestModel->first);
                                          });

            long double logDifference = log10(bestProb) - log10(secondBest->second);
            bool isUncertain = (logDifference < 5.0);

            if (isUncertain)
            {
                wordUncertain++;
                uncertainCount++;
                cout << "Sample " << sample << " - Actual: " << wordMap[actualWord].word
                     << " - Result: Uncertain (between '" << wordMap[predictedWord].word
                     << "' and '" << wordMap[secondBest->first].word
                     << "'), Log diff: " << logDifference << "\n";
            }
            else
            {
                if (predictedWord == actualWord)
                {
                    wordCorrect++;
                    correctCount++;
                }
                cout << "Sample " << sample << " - Actual: " << wordMap[actualWord].word
                     << " - Predicted: " << wordMap[predictedWord].word << "\n";
            }
        }

        double wordAccuracy = 0.0;
        if (TEST_SET_SIZE - wordUncertain > 0)
        {
            wordAccuracy = (static_cast<double>(wordCorrect) /
                             (TEST_SET_SIZE - wordUncertain)) *
                            100.0;
        }

        cout << "\nWord '" << wordMap[actualWord].word << "' Statistics:\n"
             << "Correct: " << wordCorrect << "\n"
             << "Uncertain: " << wordUncertain << "\n"
             << "Accuracy (excluding uncertain): " << fixed << setprecision(2)
             << wordAccuracy << "%\n";
    }

    double overallAccuracy = 0.0;
    if (totalTests - uncertainCount > 0)
    {
        overallAccuracy = (static_cast<double>(correctCount) /
                           (totalTests - uncertainCount)) *
                          100.0;
    }

    cout << "\nOverall Recognition Statistics:\n"
         << "Total Tests: " << totalTests << "\n"
         << "Total Correct: " << correctCount << "\n"
         << "Total Uncertain: " << uncertainCount << "\n"
         << "Overall Accuracy (excluding uncertain): " << fixed << setprecision(2)
         << overallAccuracy << "%\n";
}

void evaluateLiveRecording(const vector<HMMParameters> &wordModels)
{
    LiveRecording recorder;
    cout << "\nStarting live recording...\n";

    recorder.StartRecording();
    recorder.SaveToFile();
    recorder.ProcessRecording();
    recorder.GenerateObservationSequence();

    cout << "Processing completed. Evaluating recording...\n";

    // Get probabilities for all models with debug information
    vector<pair<int, long double>> modelProbabilities;
    cout << "\nCalculated probabilities for each word:\n";
    cout << "----------------------------------------\n";
    
    for (int modelWord = 0; modelWord < WORD_COUNT; modelWord++)
    {
        long double probability = evaluateSequence(wordModels[modelWord],
                                                   recorder.GetObservationFilePath());
        modelProbabilities.push_back(make_pair(modelWord, probability));
        
        // Print debug information
        cout << "Word '" << wordMap[modelWord].word << "': " 
             << scientific << setprecision(6) << probability << "\n";
    }

    // Check if all probabilities are valid
    bool allNegativeInf = true;
    for (size_t i = 0; i < modelProbabilities.size(); i++) {
        if (modelProbabilities[i].second > -numeric_limits<long double>::infinity()) {
            allNegativeInf = false;
            break;
        }
    }

    if (allNegativeInf) {
        cout << "\nWarning: All probabilities are extremely low. Check input audio quality.\n";
        return;
    }

    // Find best match
    vector<pair<int, long double>>::iterator bestModel = 
        max_element(modelProbabilities.begin(), modelProbabilities.end(),
                   [](const pair<int, long double> &a, const pair<int, long double> &b)
                   { return a.second < b.second; });

    // Find second best for uncertainty check
    vector<pair<int, long double>>::iterator secondBest = 
        max_element(modelProbabilities.begin(), modelProbabilities.end(),
                   [bestModel](const pair<int, long double> &a, const pair<int, long double> &b)
                   {
                       return (a.first == bestModel->first || a.second < b.second) &&
                              (b.first != bestModel->first);
                   });

    // Calculate normalized confidence scores
    long double maxLogProb = -numeric_limits<long double>::infinity();
    long double minLogProb = numeric_limits<long double>::infinity();

    // Find max and min log probabilities
    for (size_t i = 0; i < modelProbabilities.size(); i++) {
        if (modelProbabilities[i].second > maxLogProb) 
            maxLogProb = modelProbabilities[i].second;
        if (modelProbabilities[i].second < minLogProb) 
            minLogProb = modelProbabilities[i].second;
    }

    // Only calculate log difference if probabilities are valid
    long double logDifference = 0.0;
    if (bestModel->second > -numeric_limits<long double>::infinity() && 
        secondBest->second > -numeric_limits<long double>::infinity()) {
        logDifference = log10(bestModel->second) - log10(secondBest->second);
    }

    bool isUncertain = (logDifference < 5.0) || (bestModel->second < -1e10);

    cout << "\nRecognition Result:\n";
    cout << "----------------------------------------\n";
    if (isUncertain)
    {
        cout << "Result: Uncertain between '" << wordMap[bestModel->first].word 
             << "' and '" << wordMap[secondBest->first].word << "'\n"
             << "Log difference: " << logDifference << "\n"
             << "Confidence too low - please speak clearly and try again\n";
    }
    else
    {
        double confidence = (bestModel->second > maxLogProb - 10.0) ? 100.0 : 
                          exp((bestModel->second - maxLogProb)) * 100.0;
        
        cout << "Detected Command: '" << wordMap[bestModel->first].word << "'\n"
             << "Confidence: " << fixed << setprecision(2) << confidence << "%\n";
    }
}

// Main function
int main()
{
    vector<HMMParameters> wordModels(WORD_COUNT);
    bool isModelTrained = false;

    while (true)
    {
        cout << "\nHidden Markov Model Word Recognition System\n"
             << "----------------------------------------\n"
             << "1. Train New Models\n"
             << "2. Load Existing Models\n"
             << "3. Test Models\n"
             << "4. Live Recognition\n"
             << "5. Exit\n"
             << "Enter choice (1-5): ";

        char choice;
        cin >> choice;

        switch (choice)
        {
        case '1':
            cout << "\nStarting training process...\n";
            trainModels(wordModels);
            isModelTrained = true;
            break;

        case '2':
            cout << "\nLoading saved models...\n";
            isModelTrained = true;
            for (int word = 0; word < WORD_COUNT; word++)
            {
                if (!loadModelFromFile(wordModels[word], word))
                {
                    isModelTrained = false;
                    cout << "Error loading models. Please train new models.\n";
                    break;
                }
            }
            if (isModelTrained)
            {
                cout << "Models loaded successfully.\n";
            }
            break;

        case '3':
            if (!isModelTrained)
            {
                cout << "\nError: Models must be trained or loaded first.\n";
            }
            else
            {
                testModels(wordModels);
            }
            break;

        case '4':
            if (!isModelTrained)
            {
                cout << "\nError: Models must be trained or loaded first.\n";
            }
            else
            {
                evaluateLiveRecording(wordModels);
            }
            break;

        case '5':
            cout << "\nExiting program...\n";
            return 0;

        default:
            cout << "\nInvalid choice. Please try again.\n";
        }
    }
    return 0;
}