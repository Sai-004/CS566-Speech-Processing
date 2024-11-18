// 210101070_HMM1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>

using namespace std;

// Load matrix from file
void loadMatrixFromFile(const string& filepath, vector<vector<double>>& mat, int numRows, int numCols) 
{
    ifstream inputFile(filepath.c_str());
    if (!inputFile) {
        cerr << "Error: Unable to open the file " << filepath.c_str() << endl;
        exit(EXIT_FAILURE); // Exit on error
    }

    mat.resize(numRows, vector<double>(numCols));
    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            inputFile >> mat[i][j];
        }
    }
    inputFile.close();
}

// Load vector from file
void loadVectorFromFile(const string& filepath, vector<double>& vec, int length) 
{
    ifstream inputFile(filepath.c_str());
    if (!inputFile) {
        cerr << "Error: Unable to open the file " << filepath.c_str() << endl;
        exit(EXIT_FAILURE);
    }

    vec.resize(length);
    for (int i = 0; i < length; ++i) {
        inputFile >> vec[i];
    }
    inputFile.close();
}

// For printing matrix to debugging
void displayMatrix(const vector<vector<double>>& mat, const string& matrixName) 
{
    cout << "Displaying Matrix: " << matrixName.c_str() << endl;
    for (size_t i = 0; i < mat.size(); i++) {
        for (size_t j = 0; j < mat[i].size(); j++) {
            cout << mat[i][j] << " \t";
        }
        cout << endl;
    }
	cout << endl;
}

// Forward Algorithm implementation
double computeForwardProbability(const vector<vector<double>>& A, const vector<vector<double>>& B, const vector<double>& init_Pi, const vector<int>& ObservationSequence) 
{
    int numStates = A.size(); // Number of states (N)
    int sequenceLength = ObservationSequence.size(); // Observation sequence length (T)

    vector<vector<double>> AlphaMatrix(sequenceLength, vector<double>(numStates, 0.0));

    // Initialization process
    for (int i = 0; i < numStates; i++) {
        AlphaMatrix[0][i] = init_Pi[i] * B[i][ObservationSequence[0]];
    }

    // Induction process for Alpha 
    for (int t = 1; t < sequenceLength; t++) {
        for (int j = 0; j < numStates; j++) {
            double s = 0.0;
            for (int i = 0; i < numStates; i++) {
                s += AlphaMatrix[t-1][i] * A[i][j];
            }
            AlphaMatrix[t][j] = s * B[j][ObservationSequence[t]];
        }
    }

    displayMatrix(AlphaMatrix, "Alpha Matrix (Forward Probability Matrix)");

    // Termination process
    double finalProb = 0.0;
    for (int i = 0; i < numStates; ++i) {
        finalProb += AlphaMatrix[sequenceLength-1][i];
    }

    return finalProb;
}

// Backward Algorithm implementation
double computeBackwardProbability(const vector<vector<double>>& A, const vector<vector<double>>& B, const vector<double>& init_Pi, const vector<int>& ObservationSequence) 
{
    int numStates = A.size();  // Number of states
    int sequenceLength = ObservationSequence.size();  // Length of the observation sequence

    vector<vector<double>> BetaMatrix(sequenceLength, vector<double>(numStates, 0.0));

    // Initializing process
    for (int i = 0; i < numStates; i++) {
        BetaMatrix[sequenceLength - 1][i] = 1.0;
    }

    // Induction process for Beta
    for (int t = sequenceLength - 2; t >= 0; t--) {
        for (int i = 0; i < numStates; ++i) {
            double s = 0.0;
            for (int j = 0; j < numStates; j++) {
                s += A[i][j] * B[j][ObservationSequence[t + 1]] * BetaMatrix[t + 1][j];
            }
            BetaMatrix[t][i] = s;
        }
    }

    displayMatrix(BetaMatrix, "Beta Matrix (Backward Probability Matrix)");

    // Termination process
    double finalProb = 0.0;
    for (int i = 0; i < numStates; ++i) {
        finalProb += init_Pi[i] * B[i][ObservationSequence[0]] * BetaMatrix[0][i];
    }

    return finalProb;
}

int main() {
    vector<vector<double>> transitionMatrix; // A: State transition probabilities
    vector<vector<double>> observationMatrix; // B: Observation probabilities
    vector<double> initialProbabilities;  // Pi: Initial state distribution
    vector<int> observationSequence;  // O: Sequence of ObservationSequence

    // Load matrices and vectors from files
    loadMatrixFromFile("files/A_ijMatrix.txt", transitionMatrix, 5, 5);
    loadMatrixFromFile("files/B_jkMatrix.txt", observationMatrix, 5, 32);
    loadVectorFromFile("files/P_iVector.txt", initialProbabilities, 5);

    // Read observation sequence from file
    ifstream observationFile("files/ObservationSequence.txt");
    if (!observationFile) {
        cerr << "Error: Unable to open observation sequence file" << endl;
        return EXIT_FAILURE;
    }

    int observationValue;
    while (observationFile >> observationValue) {
        observationSequence.push_back(observationValue - 1);  // Convert to zero-based index
    }
    observationFile.close();

    // Calculate forward and backward probabilities
    double forwardProb = computeForwardProbability(transitionMatrix, observationMatrix, initialProbabilities, observationSequence);
    double backwardProb = computeBackwardProbability(transitionMatrix, observationMatrix, initialProbabilities, observationSequence);

    // Output the results
    cout << "Forward Probability of the Observation Sequence: " << forwardProb << endl;
    cout << "Backward Probability of the Observation Sequence: " << backwardProb << endl;

    return 0;
}