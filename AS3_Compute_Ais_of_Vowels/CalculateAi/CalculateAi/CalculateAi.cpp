// CalculateAi.cpp : Defines the entry point for the console application.

#include "stdafx.h"
#include <iostream>
#include <math.h>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <sys/stat.h>

using namespace std;

#define FRAME_SIZE 320 // Assuming 320 samples/frame
#define p 12  // Order of LPC
#define PI 3.14159265358979323846 // Value of pi for calculation
#define MAX_SIZE 10000 // Max data size

// Computing autocorrelation coefficients R_i
void ComputeRi(double* frame, int frameSize, double* Ri) {
    for (int i = 0; i <= p; ++i) {
        Ri[i] = 0.0;
        for (int j = 0; j < frameSize - i; ++j) {
            Ri[i] += frame[j] * frame[j + i];
        }
    }
}

// Computing LPC coefficients using Levinson-Durbin algorithm 
void ComputeAi(double* R, double* A) {
    double alpha[p + 1][p + 1] = {0};  // 2D array for LPC coefficients
    double E[p + 1] = {0};  // Error predictions
    double K[p + 1] = {0};  // Reflection coefficients

    E[0] = R[0];

    for (int i = 1; i <= p; ++i) {
        double sum = 0.0;
        for (int j = 1; j < i; ++j) {
            sum += alpha[i - 1][j] * R[i - j];
        }

        K[i] = (R[i] - sum) / E[i - 1];
        alpha[i][i] = K[i];  // alpha(i)_i = K_i

        for (int j = 1; j < i; ++j) {
            alpha[i][j] = alpha[i - 1][j] - K[i] * alpha[i - 1][i - j];
        }

        E[i] = (1.0 - K[i] * K[i]) * E[i - 1];
    }

    // LPC coefficients (last row in the 2D array)
    for (int j = 1; j <= p; ++j) {
        A[j] = alpha[p][j];
    }
}

// Performing DC shift and Amplitude normalizations
void applyNormalization(double* data, int dataSize) {
    // DC shift normalization
    double totalShift = 0.0;
    for (int i = 0; i < dataSize; ++i) {
        totalShift += data[i];
    }
    double DCShift = totalShift / dataSize;

    for (int i = 0; i < dataSize; ++i) {
        data[i] -= DCShift;
    }

    // Amplitude normalization
    double maxAmplitude = 0.0;
    for (int i = 0; i < dataSize; ++i) {
        if (fabs(data[i]) > maxAmplitude) {
            maxAmplitude = fabs(data[i]);
        }
    }

    if (maxAmplitude != 0) {
        for (int i = 0; i < dataSize; ++i) {
            data[i] /= maxAmplitude;
        }
    }
}

// Applying Hamming Window to data frame
void applyHammingWindow(double* frame, int frameSize) {
    for (int i = 0; i < frameSize; ++i) {
        frame[i] *= 0.54 - 0.46 * cos(2 * PI * i / (frameSize - 1));
    }
}

// Calculating energy of a frame
double calculateFrameEnergy(double* frame, int frameSize) {
    double energy = 0.0;
    for (int i = 0; i < frameSize; ++i) {
        energy += frame[i] * frame[i];
    }
    return energy;
}

// Finding frame with maximum energy
int findMaxEnergyFrame(vector<double>& frameEnergies) {
    return max_element(frameEnergies.begin(), frameEnergies.end()) - frameEnergies.begin();
}

void processFile(const string& filename) {
	// Reading the data from the input file
    ifstream inputFile("Samples/" + filename);
    if (!inputFile.is_open()) {
        cerr << "Error: Could not open the input file: " << filename << endl;
        return;
    }

    vector<double> amplitudes;
    double value;
    while (inputFile >> value) {
        amplitudes.push_back(value);
    }
    inputFile.close();

    int dataSize = amplitudes.size();

    // Applying normalizations on entire data
    applyNormalization(amplitudes.data(), dataSize);

    // Prepare output file paths
    string aiOutputPath = "Results/Ai_values/" + filename + "_ai.txt";
    string riOutputPath = "Results/Ri_values/" + filename + "_ri.txt";

    ofstream outputA(aiOutputPath);
    ofstream outputR(riOutputPath);

    if (!outputA.is_open() || !outputR.is_open()) {
        cerr << "Error: Could not open the output files for " << filename << endl;
        return;
    }

    // Calculate energies for all frames
    vector<double> frameEnergies;
    for (int i = 0; i + FRAME_SIZE <= dataSize; i += FRAME_SIZE) {
        double energy = calculateFrameEnergy(&amplitudes[i], FRAME_SIZE);
        frameEnergies.push_back(energy);
    }

    // Find the frame with maximum energy
    int maxEnergyFrame = findMaxEnergyFrame(frameEnergies);

    // Process 5 frames around the max energy frame
    for (int frameOffset = -2; frameOffset <= 2; ++frameOffset) {
        int currentFrame = maxEnergyFrame + frameOffset;
        
        // Ensure we're not going out of bounds
        if (currentFrame < 0 || currentFrame >= frameEnergies.size()) {
            continue;
        }

        int startIndex = currentFrame * FRAME_SIZE;
        vector<double> frame(amplitudes.begin() + startIndex, amplitudes.begin() + startIndex + FRAME_SIZE);

        // Applying Hamming Window to the frame
        applyHammingWindow(frame.data(), FRAME_SIZE);

        // Computing autocorrelation coefficients R_i
        double Ri[p + 1];
        ComputeRi(frame.data(), FRAME_SIZE, Ri);

        // R coefficients
        outputR << "Frame " << currentFrame << " (offset " << frameOffset << "): ";
        for (int j = 0; j <= p; ++j) {
            outputR << Ri[j] << "\t";
        }
        outputR << endl;

        // Computing LPC coeffs alpha(p)_j using Levinson Durbin Algorithm
        double A[p + 1];
        ComputeAi(Ri, A);

        // A coefficients
        outputA << "Frame " << currentFrame << " (offset " << frameOffset << "): ";
        for (int j = 1; j <= p; ++j) {
            outputA << A[j] << "\t";
        }
        outputA << endl;
    }

    outputA.close();
    outputR.close();
    cout << "Processed file: " << filename << endl;
}

// Main function: Processing data and computing R and A coefficients
int main() {
	vector<string> vowels;
    vowels.push_back("a");
    vowels.push_back("e");
    vowels.push_back("i");
    vowels.push_back("o");
    vowels.push_back("u");

    // Process each vowel and each file
    for (size_t v = 0; v < vowels.size(); ++v) {
        string vowel = vowels[v];
        for (int i = 1; i <= 20; ++i) {
            // Generate the input filename
            stringstream filenameStream;
            filenameStream << "210101070_" << vowel << "_" << i << ".txt";
            string inputFilename = filenameStream.str();

            // Process the file
            processFile(inputFilename);
        }
    }

	cout << "Processing complete. Check 'Results' folder for A and R coefficients of all processed files." << endl;
    
    return 0;
}