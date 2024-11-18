#include <iostream>
#include <string>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <numeric>

using namespace std;
namespace fs = filesystem;

const double PI = 3.14159265358979323846;
const int SEGMENT_WIDTH = 320;
const int SEGMENT_STRIDE = 240;
const int PREDICTION_ORDER = 12;
const int MAX_SIGNAL_LENGTH = 150000;

class SignalProcessor
{
private:
    // Perform signal baseline adjustment
    void centerSignal(vector<double> &signal)
    {
        double totalSum = accumulate(signal.begin(), signal.end(), 0.0);
        double signalMean = totalSum / signal.size();

        double maxAmplitude = 0.0;
        for (auto &sample : signal)
        {
            sample -= signalMean;
            maxAmplitude = max(maxAmplitude, abs(sample));
        }

        if (maxAmplitude > 0)
        {
            for (auto &sample : signal)
            {
                sample /= maxAmplitude;
            }
        }
    }

    // Apply windowing function
    void applySpectralWindow(vector<double> &segment)
    {
        for (size_t i = 0; i < segment.size(); ++i)
        {
            segment[i] *= 0.54 - 0.46 * cos(2 * PI * i / (segment.size() - 1));
        }
    }

    // Compute signal cross-correlation
    vector<double> calculateCorrelation(const vector<double> &segment, int predOrder)
    {
        vector<double> correlations(predOrder + 1, 0.0);
        for (int i = 0; i <= predOrder; ++i)
        {
            for (size_t j = 0; j < segment.size() - i; ++j)
            {
                correlations[i] += segment[j] * segment[j + i];
            }
        }
        return correlations;
    }

    // Advanced linear prediction coefficients estimation
    vector<double> estimatePredictionCoefficients(const vector<double> &correlations)
    {
        vector<double> predCoeffs(PREDICTION_ORDER + 1, 0.0);
        vector<vector<double>> alphaMatrix(PREDICTION_ORDER + 1, vector<double>(PREDICTION_ORDER + 1, 0.0));
        vector<double> errorSignal(PREDICTION_ORDER + 1, 0.0);
        vector<double> reflectionCoeffs(PREDICTION_ORDER + 1, 0.0);

        errorSignal[0] = correlations[0];

        for (int i = 1; i <= PREDICTION_ORDER; ++i)
        {
            double interimSum = 0.0;
            for (int j = 1; j < i; ++j)
            {
                interimSum += alphaMatrix[i - 1][j] * correlations[i - j];
            }

            reflectionCoeffs[i] = (correlations[i] - interimSum) / errorSignal[i - 1];
            alphaMatrix[i][i] = reflectionCoeffs[i];

            for (int j = 1; j < i; ++j)
            {
                alphaMatrix[i][j] = alphaMatrix[i - 1][j] - reflectionCoeffs[i] * alphaMatrix[i - 1][i - j];
            }

            errorSignal[i] = (1.0 - reflectionCoeffs[i] * reflectionCoeffs[i]) * errorSignal[i - 1];
        }

        for (int j = 1; j <= PREDICTION_ORDER; ++j)
        {
            predCoeffs[j] = alphaMatrix[PREDICTION_ORDER][j];
        }

        return predCoeffs;
    }

    // Compute spectral features
    vector<double> calculateSpectralFeatures(const vector<double> &predCoeffs)
    {
        vector<double> spectralFeats(PREDICTION_ORDER + 1, 0.0);

        for (int m = 1; m <= PREDICTION_ORDER; ++m)
        {
            spectralFeats[m] = predCoeffs[m];
            for (int k = 1; k < m; ++k)
            {
                spectralFeats[m] += (k * spectralFeats[k] * predCoeffs[m - k]) / m;
            }
        }

        return spectralFeats;
    }

    // Apply spectral enhancement window
    void enhanceSpectralFeatures(vector<double> &spectralFeats)
    {
        for (size_t i = 1; i < spectralFeats.size(); ++i)
        {
            spectralFeats[i] *= (1.0 + (PREDICTION_ORDER / 2.0) * sin(PI * i / PREDICTION_ORDER));
        }
    }

public:
    struct SignalSegment
    {
        vector<double> rawData;
        vector<double> correlationValues;
        vector<double> predictionCoefficients;
        vector<double> spectralFeatures;
    };

    vector<SignalSegment> processSignal(const string &filePath)
    {
        vector<SignalSegment> signalSegments;
        vector<double> rawSignalData = readSignalFromFile(filePath);

        centerSignal(rawSignalData);

        // Segment the signal
        for (size_t i = 0; i + SEGMENT_WIDTH <= rawSignalData.size(); i += SEGMENT_STRIDE)
        {
            SignalSegment segment;
            segment.rawData.assign(rawSignalData.begin() + i, rawSignalData.begin() + i + SEGMENT_WIDTH);

            applySpectralWindow(segment.rawData);
            segment.correlationValues = calculateCorrelation(segment.rawData, PREDICTION_ORDER);
            segment.predictionCoefficients = estimatePredictionCoefficients(segment.correlationValues);
            segment.spectralFeatures = calculateSpectralFeatures(segment.predictionCoefficients);

            enhanceSpectralFeatures(segment.spectralFeatures);

            signalSegments.push_back(segment);
        }

        return signalSegments;
    }

    // Utility function to read signal from file
    vector<double> readSignalFromFile(const string &filepath)
    {
        vector<double> signalData;
        ifstream inputFile(filepath);
        double value;

        while (inputFile >> value)
        {
            signalData.push_back(value);
        }

        return signalData;
    }

    // Write results to CSV
    void exportResultsToCSV(const string &outputFilename,
                            const vector<string> &fileList,
                            const vector<vector<SignalSegment>> &processedSignals)
    {
        ofstream csvFile(outputFilename);

        // Write header
        csvFile << "Filename,SegmentIndex";
        for (int i = 1; i <= PREDICTION_ORDER; ++i)
        {
            csvFile << ",Feature" << i;
        }
        csvFile << "\n";

        // Write data
        for (size_t i = 0; i < fileList.size(); ++i)
        {
            for (size_t j = 0; j < processedSignals[i].size(); ++j)
            {
                csvFile << fileList[i] << "," << j;
                for (size_t k = 1; k <= PREDICTION_ORDER; ++k)
                {
                    csvFile << "," << processedSignals[i][j].spectralFeatures[k];
                }
                csvFile << "\n";
            }
        }

        csvFile.close();
    }

    // Discover text files in a directory
    vector<string> discoverTextFiles(const string &directoryPath)
    {
        vector<string> fileList;
        for (const auto &entry : fs::directory_iterator(directoryPath))
        {
            if (entry.is_regular_file() && entry.path().extension() == ".txt")
            {
                fileList.push_back(entry.path().filename().string());
            }
        }
        return fileList;
    }
};

int main()
{
    SignalProcessor processor;
    string inputDirectory = "210101070_dataset/English/txt";
    string outputFile = "CepstralCoefficients.csv";

    cout << "Processing files from: " << inputDirectory << endl;

    auto fileList = processor.discoverTextFiles(inputDirectory);
    vector<vector<SignalProcessor::SignalSegment>> processedSignals;

    for (size_t i = 0; i < fileList.size(); ++i)
    {
        cout << "Processing file " << (i + 1) << " of " << fileList.size() << ": " << fileList[i] << endl;
        string fullPath = inputDirectory + "/" + fileList[i];
        auto segmentedSignals = processor.processSignal(fullPath);
        processedSignals.push_back(segmentedSignals);
    }

    processor.exportResultsToCSV(outputFile, fileList, processedSignals);

    cout << "Processing complete. Results written to: " << outputFile << endl;

    return 0;
}