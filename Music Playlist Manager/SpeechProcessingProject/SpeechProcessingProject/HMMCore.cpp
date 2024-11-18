#include "stdafx.h"
#include "HMMCore.h"
#include <fstream>
#include <string>
#include <limits>

struct ComputationMatrices {
    int observationSeq[MAX_SEQUENCE_LENGTH];
    int seqLength;
    long double forwardProbs[MAX_SEQUENCE_LENGTH][STATE_COUNT];
};

static ComputationMatrices compute;

bool loadModelFromFile(HMMParameters& model, const std::string& filepath) {
    std::ifstream file(filepath.c_str());
    if (!file) return false;

    std::string line;
    std::getline(file, line); // Skip header

    for (int i = 0; i < STATE_COUNT; i++) {
        for (int j = 0; j < STATE_COUNT; j++) {
            file >> model.transitionMatrix[i][j];
        }
    }

    std::getline(file, line); // Skip newline
    std::getline(file, line); // Skip header
    std::getline(file, line); // Skip newline

    for (int i = 0; i < STATE_COUNT; i++) {
        for (int j = 0; j < OBSERVATION_SYMBOLS; j++) {
            file >> model.emissionMatrix[i][j];
        }
    }

    std::getline(file, line); // Skip newline
    std::getline(file, line); // Skip header
    std::getline(file, line); // Skip newline

    for (int i = 0; i < STATE_COUNT; i++) {
        file >> model.initialProbs[i];
    }

    file.close();
    return true;
}

bool loadObservationSequence(const std::string& filepath) {
    std::ifstream file(filepath.c_str());
    if (!file) return false;

    compute.seqLength = 0;
    int observation;
    while (file >> observation && compute.seqLength < MAX_SEQUENCE_LENGTH) {
        compute.observationSeq[compute.seqLength++] = observation;
    }
    file.close();
    return true;
}

long double calculateSequenceProbability(const HMMParameters& model) {
    // Initialize first time step
    for (int i = 0; i < STATE_COUNT; i++) {
        compute.forwardProbs[0][i] = model.initialProbs[i] * 
                                    model.emissionMatrix[i][compute.observationSeq[0]];
    }

    // Forward recursion
    for (int t = 1; t < compute.seqLength; t++) {
        for (int j = 0; j < STATE_COUNT; j++) {
            compute.forwardProbs[t][j] = 0.0;
            for (int i = 0; i < STATE_COUNT; i++) {
                compute.forwardProbs[t][j] += compute.forwardProbs[t-1][i] * 
                                            model.transitionMatrix[i][j];
            }
            compute.forwardProbs[t][j] *= model.emissionMatrix[j][compute.observationSeq[t]];
        }
    }

    // Calculate total probability
    long double probability = 0.0;
    for (int i = 0; i < STATE_COUNT; i++) {
        probability += compute.forwardProbs[compute.seqLength-1][i];
    }
    return probability;
}

long double evaluateSequence(const HMMParameters& model, const std::string& filepath) {
    if (!loadObservationSequence(filepath)) {
        return -std::numeric_limits<long double>::infinity();
    }
    return calculateSequenceProbability(model);
}