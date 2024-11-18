#pragma once
#include "HMMCore.h"
#include <string>
#include <vector>

// Forward declare Form1
namespace SpeechProcessingProject {
    ref class Form1;
}

// Define the Test class
class Test {
public:
    static bool GenerateTestObservations(const std::string& testDataDir);
    static bool RunCommandRecognitionTest(HMMParameters* commandModels, const std::string& testDataDir);
    static bool RunDigitRecognitionTest(HMMParameters* digitModels, const std::string& testDataDir);
    static bool RunCompleteFlowTest(SpeechProcessingProject::Form1^ form, const std::string& testDataDir);
    
private:
    static bool CreateTestFile(const std::string& filename, const std::vector<int>& observations);
    static std::vector<int> LoadObservations(const std::string& filename);

    // Static arrays for test data
    static const int commandSamples[][100];
    static const int digitSamples[][100];
    static const int commandSampleLengths[];
    static const int digitSampleLengths[];
};