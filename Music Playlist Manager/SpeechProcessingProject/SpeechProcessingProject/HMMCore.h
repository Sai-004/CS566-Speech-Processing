#pragma once

// Configuration constants
const int STATE_COUNT = 5;
const int MAX_SEQUENCE_LENGTH = 100;
const int OBSERVATION_SYMBOLS = 32;
const double MIN_PROBABILITY = 1e-15;

// HMM model structure
struct HMMParameters {
    double transitionMatrix[STATE_COUNT][STATE_COUNT];
    double emissionMatrix[STATE_COUNT][OBSERVATION_SYMBOLS];
    double initialProbs[STATE_COUNT];
};

// Function declarations
bool loadModelFromFile(HMMParameters& model, const std::string& filepath);
long double evaluateSequence(const HMMParameters& model, const std::string& filepath);