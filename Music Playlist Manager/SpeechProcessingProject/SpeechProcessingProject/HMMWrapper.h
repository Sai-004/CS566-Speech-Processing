// HMMWrapper.h
#pragma once
#include "HMMCore.h"
#include "LiveRecording.h"
#include <string>
#include <vector>

class HMMWrapper {
private:
    LiveRecording* recorder;
    std::vector<HMMParameters> commandModels;
    std::vector<HMMParameters> digitModels;
    const char* commands[5] = {"create", "add", "play", "stop", "shuffle"};

public:
    HMMWrapper() {
        recorder = new LiveRecording();
        commandModels.resize(5);
        digitModels.resize(10);
    }

    ~HMMWrapper() {
        if (recorder) {
            delete recorder;
            recorder = nullptr;
        }
    }

    bool LoadModels(const std::string& modelPath) {
        try {
            // Load command models
            for (int i = 0; i < 5; i++) {
                std::string modelFile = modelPath + "\\command_" + commands[i] + ".txt";
                if (!loadModelFromFile(commandModels[i], modelFile))
                    return false;
            }

            // Load digit models
            for (int i = 0; i < 10; i++) {
                char digit[2] = { '0' + i, 0 };
                std::string modelFile = modelPath + "\\digit_" + digit + ".txt";
                if (!loadModelFromFile(digitModels[i], modelFile))
                    return false;
            }
            return true;
        }
        catch (...) {
            return false;
        }
    }

    void StartRecording() {
        if (recorder) recorder->StartRecording();
    }

    void ProcessRecording() {
        if (recorder) {
            recorder->SaveToFile();
            recorder->ProcessRecording();
            recorder->GenerateObservationSequence();
        }
    }

    const char* GetObservationPath() {
        return recorder ? recorder->GetObservationFilePath() : "";
    }

    int RecognizeCommand() {
        if (!recorder) return -1;

        std::string obsPath = recorder->GetObservationFilePath();
        std::vector<std::pair<int, long double>> modelProbabilities;

        for (int i = 0; i < 5; i++) {
            long double probability = evaluateSequence(commandModels[i], obsPath);
            modelProbabilities.push_back(std::make_pair(i, probability));
        }

        auto bestMatch = std::max_element(modelProbabilities.begin(), 
            modelProbabilities.end(),
            [](const std::pair<int, long double>& a, 
               const std::pair<int, long double>& b) {
                return a.second < b.second;
            });

        return bestMatch->first;
    }

    int RecognizeDigit() {
        if (!recorder) return -1;

        std::string obsPath = recorder->GetObservationFilePath();
        std::vector<std::pair<int, long double>> modelProbabilities;

        for (int i = 0; i < 10; i++) {
            long double probability = evaluateSequence(digitModels[i], obsPath);
            modelProbabilities.push_back(std::make_pair(i, probability));
        }

        auto bestMatch = std::max_element(modelProbabilities.begin(), 
            modelProbabilities.end(),
            [](const std::pair<int, long double>& a, 
               const std::pair<int, long double>& b) {
                return a.second < b.second;
            });

        return bestMatch->first;
    }

    const char* GetCommandName(int index) {
        if (index >= 0 && index < 5)
            return commands[index];
        return "";
    }
};