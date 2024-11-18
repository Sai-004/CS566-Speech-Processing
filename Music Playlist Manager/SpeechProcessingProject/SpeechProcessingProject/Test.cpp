#include "stdafx.h"
#include "Test.h"
#include "Form1.h"
#include <iostream>
#include <fstream>
#include <string>

// Static arrays for command samples
const int Test::commandSampleLengths[] = {75, 79, 75, 79, 79};
const int Test::commandSamples[][100] = {
    // create (75)
    {1, 6, 6, 1, 6, 7, 4, 7, 6, 1, 6, 6, 1, 6, 6, 1, 1, 6, 4, 6, 1, 1, 7, 4, 7, 26, 26, 14, 11, 1, 15, 8, 8, 8, 11, 11, 8, 8, 11, 11, 8, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 10, 10, 10, 10, 10, 29, 29, 29, 1, 1, 5, 1, 1, 1, 6, 1, 1, 1, 4, 7, 7, 1, 1},
    // add (79)
    {0, 0, 0, 2, 2, 3, 3, 3, 3, 4, 0, 0, 0, 0, 0, 0, 0, 3, 28, 28, 28, 28, 28, 28, 28, 28, 3, 25, 24, 25, 25, 25, 25, 25, 24, 25, 25, 25, 25, 25, 25, 15, 15, 15, 15, 14, 14, 14, 14, 14, 11, 11, 11, 14, 10, 10, 1, 1, 1, 1, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 23, 23, 23, 23, 0, 4, 0},
    // play (75)
    {2, 2, 3, 2, 1, 20, 16, 20, 20, 20, 2, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 11, 11, 9, 11, 11, 11, 2, 29, 17, 29, 29, 29, 29, 29, 29, 29, 1, 6, 4, 4, 0, 1, 4, 1, 1, 1, 1, 1},
    // stop (79)
    {3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 3, 7, 6, 22, 22, 22, 22, 20, 21, 20, 20, 20, 20, 20, 20, 21, 22, 16, 22, 20, 20, 20, 20, 19, 17, 17, 17, 17, 17, 17, 17, 16, 16, 16, 16, 16, 16, 17, 28, 28, 29, 29, 29, 29, 31, 5, 6, 5, 7, 7, 5, 7, 5, 6, 6, 6, 6, 6, 6, 5, 6, 5, 6, 6, 5, 5, 3, 29, 29},
    // shuffle (79)
    {6, 6, 6, 7, 5, 5, 6, 6, 6, 6, 6, 5, 1, 5, 6, 6, 6, 6, 6, 6, 7, 18, 18, 31, 31, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 11, 30, 14, 30, 30, 19, 30, 30, 30, 30, 30, 19, 30, 15, 15, 30, 30, 23, 23, 23, 23, 23, 29, 29, 3, 28, 18, 18, 18, 18, 18, 18, 18, 18, 31, 31, 31, 6, 5, 5, 27, 15, 15, 15}
};

// Static arrays for digit samples
const int Test::digitSampleLengths[] = {79, 75, 79, 79, 79, 79, 79, 79, 79, 79};
const int Test::digitSamples[][100] = {
    // 0
    {28, 30, 30, 30, 30, 31, 30, 7, 23, 22, 31, 31, 30, 23, 31, 31, 31, 20, 21, 7, 21, 21, 20, 31, 7, 5, 5, 9, 9, 9, 9, 9, 9, 17, 17, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 3, 15, 13, 13, 13, 13, 13, 13, 13, 13, 12, 12, 13, 13, 13, 11, 11, 11, 11, 11, 11, 6, 20, 30, 31, 7, 21, 21, 20, 20, 31, 7, 31, 31},
    
    // 1
    {18, 18, 20, 20, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 4, 4, 10, 10, 8, 8, 12, 12, 14, 14, 14, 14, 14, 23, 15, 15, 15, 3, 15, 3, 2, 2, 2, 2, 18, 2, 2, 2, 2, 18, 18, 18, 2, 18, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6},
    
    // 2
    {21, 23, 23, 29, 29, 29, 29, 29, 23, 30, 23, 20, 22, 20, 20, 22, 7, 21, 21, 20, 20, 20, 7, 20, 20, 21, 20, 21, 20, 22, 20, 20, 22, 7, 21, 21, 20, 31, 19, 15, 13, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 6, 20, 20, 22, 21, 21, 11, 4, 6, 20, 11, 4, 10, 4, 4, 4},
    
    // 3
    {21, 20, 31, 31, 7, 31, 31, 21, 20, 21, 31, 22, 31, 20, 31, 7, 21, 28, 28, 28, 30, 7, 31, 31, 21, 20, 21, 31, 22, 31, 31, 13, 27, 17, 15, 15, 3, 3, 18, 3, 3, 3, 1, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 2, 6, 30, 22, 31, 31, 22, 21, 21, 21, 7, 31, 31, 7, 7, 31, 21, 31, 23, 31},
    
    // 4
    {20, 7, 20, 20, 23, 7, 21, 7, 20, 20, 7, 7, 22, 7, 20, 20, 22, 20, 7, 20, 20, 21, 7, 21, 21, 21, 23, 28, 28, 28, 29, 28, 28, 19, 10, 12, 12, 12, 8, 12, 12, 12, 13, 12, 13, 13, 19, 13, 15, 15, 19, 18, 19, 19, 18, 18, 18, 11, 11, 11, 11, 6, 6, 6, 6, 20, 20, 21, 20, 22, 20, 20, 20, 7, 21, 21, 20, 20, 20},
    
    // 5
    {31, 21, 20, 21, 20, 22, 20, 20, 20, 7, 21, 21, 20, 20, 20, 7, 20, 20, 23, 29, 29, 29, 29, 14, 14, 29, 24, 29, 29, 23, 29, 29, 23, 11, 13, 13, 13, 13, 13, 15, 15, 15, 3, 3, 3, 3, 2, 2, 0, 2, 0, 2, 2, 2, 2, 10, 10, 10, 10, 11, 30, 21, 21, 7, 20, 20, 7, 7, 20, 21, 20, 21, 20, 22, 20, 20, 22, 21, 21},
    
    // 6
    {20, 31, 7, 7, 31, 21, 20, 21, 31, 22, 23, 31, 22, 21, 21, 21, 7, 20, 31, 7, 7, 31, 21, 20, 21, 31, 22, 31, 24, 24, 25, 25, 24, 25, 25, 24, 25, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 2, 30, 30, 30, 26, 24, 24, 25, 24, 24, 25, 25, 24, 24, 24, 24, 24, 24, 24, 24, 7, 20, 31, 7, 7, 31, 21, 20, 20, 22, 22, 7},
    
    // 7
    {5, 4, 4, 5, 6, 20, 30, 31, 23, 21, 21, 20, 31, 31, 7, 31, 31, 21, 31, 23, 31, 22, 31, 23, 7, 23, 21, 21, 20, 31, 31, 7, 31, 31, 21, 31, 23, 31, 22, 26, 26, 21, 26, 21, 21, 7, 31, 31, 7, 31, 31, 21, 31, 23, 31, 22, 23, 31, 22, 23, 21, 21, 7, 31, 31, 7, 7, 31, 21, 31, 23, 31, 22, 23, 31, 22, 23, 21, 21},
    
    // 8
    {2, 15, 19, 30, 15, 18, 7, 7, 28, 7, 20, 20, 7, 21, 21, 20, 20, 20, 7, 20, 20, 20, 20, 21, 20, 22, 20, 20, 20, 7, 21, 21, 20, 20, 20, 7, 20, 20, 20, 18, 16, 16, 16, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 15, 23, 15, 15, 30, 17, 15, 15, 7, 20, 20, 20, 21, 20, 21, 20, 22, 20, 20},
    
    // 9
    {17, 6, 18, 18, 18, 18, 7, 21, 21, 20, 20, 20, 7, 20, 20, 20, 20, 21, 20, 22, 20, 20, 20, 7, 21, 21, 20, 20, 20, 7, 20, 20, 21, 6, 4, 4, 6, 6, 21, 24, 24, 19, 19, 19, 3, 3, 3, 1, 3, 3, 3, 3, 2, 2, 2, 2, 2, 4, 2, 4, 4, 4, 27, 4, 5, 5, 5, 5, 6, 20, 22, 20, 20, 22, 7, 21, 21, 7, 20}
};

bool Test::CreateTestFile(const std::string& filename, const std::vector<int>& observations) {
    std::ofstream file(filename.c_str());
    if (!file) return false;

    for (size_t i = 0; i < observations.size(); i++) {
        file << observations[i];
        if (i < observations.size() - 1) file << " ";
    }
    file.close();
    return true;
}

std::vector<int> Test::LoadObservations(const std::string& filename) {
    std::vector<int> observations;
    std::ifstream file(filename.c_str());
    int value;
    while (file >> value) {
        observations.push_back(value);
    }
    return observations;
}

bool Test::GenerateTestObservations(const std::string& testDataDir) {
    static const char* commands[] = {"create", "add", "play", "stop", "shuffle"};
    
    // Create test files for commands
    for (int i = 0; i < 5; i++) {
        std::string filename = testDataDir + "\\test_command_" + commands[i] + ".txt";
        std::vector<int> observations(commandSamples[i], 
                                    commandSamples[i] + commandSampleLengths[i]);
        if (!CreateTestFile(filename, observations)) {
            return false;
        }
    }

    // Create test files for digits
    char buffer[8];
    for (int i = 0; i < 10; i++) {
        _itoa_s(i, buffer, sizeof(buffer), 10);
        std::string filename = testDataDir + "\\test_digit_" + buffer + ".txt";
        std::vector<int> observations(digitSamples[i], 
                                    digitSamples[i] + digitSampleLengths[i]);
        if (!CreateTestFile(filename, observations)) {
            return false;
        }
    }

    return true;
}

bool Test::RunCommandRecognitionTest(HMMParameters* commandModels, const std::string& testDataDir) {
    static const char* commands[] = {"create", "add", "play", "stop", "shuffle"};
    bool allTestsPassed = true;

    std::cout << "\nRunning Command Recognition Tests..." << std::endl;

    for (int i = 0; i < 5; i++) {
        std::string testFile = testDataDir + "\\test_command_" + commands[i] + ".txt";
        
        double maxProb = -std::numeric_limits<double>::infinity();
        int bestMatch = -1;

        for (int j = 0; j < 5; j++) {
            double prob = evaluateSequence(commandModels[j], testFile);
            if (prob > maxProb) {
                maxProb = prob;
                bestMatch = j;
            }
        }

        bool passed = (bestMatch == i);
        allTestsPassed &= passed;

        std::cout << "Test " << commands[i] << ": " 
                 << (passed ? "PASSED" : "FAILED")
                 << " (Recognized as: " << commands[bestMatch] << ")" << std::endl;
    }

    return allTestsPassed;
}

bool Test::RunDigitRecognitionTest(HMMParameters* digitModels, const std::string& testDataDir) {
    bool allTestsPassed = true;

    std::cout << "\nRunning Digit Recognition Tests..." << std::endl;

    char buffer[8];
    for (int i = 0; i < 10; i++) {
        _itoa_s(i, buffer, sizeof(buffer), 10);
        std::string testFile = testDataDir + "\\test_digit_" + buffer + ".txt";
        
        double maxProb = -std::numeric_limits<double>::infinity();
        int bestMatch = -1;

        for (int j = 0; j < 10; j++) {
            double prob = evaluateSequence(digitModels[j], testFile);
            if (prob > maxProb) {
                maxProb = prob;
                bestMatch = j;
            }
        }

        bool passed = (bestMatch == i);
        allTestsPassed &= passed;

        std::cout << "Test digit " << i << ": " 
                 << (passed ? "PASSED" : "FAILED")
                 << " (Recognized as: " << bestMatch << ")" << std::endl;
    }

    return allTestsPassed;
}

bool Test::RunCompleteFlowTest(SpeechProcessingProject::Form1^ form, const std::string& testDataDir) {
    bool success = true;
    std::cout << "\nRunning Complete Flow Test..." << std::endl;

    try {
        static const char* testFiles[] = {
            "test_command_create.txt",
            "test_command_add.txt",
            "test_digit_1.txt",
            "test_command_play.txt",
            "test_command_stop.txt",
            "test_command_shuffle.txt"
        };

        for (int i = 0; i < 6; i++) {
            std::string fullPath = testDataDir + "\\" + testFiles[i];
            {
                std::ifstream src(fullPath.c_str(), std::ios::binary);
                if (!src) {
                    std::cout << "Failed to open source file: " << fullPath << std::endl;
                    return false;
                }
                std::ofstream dst("live_observation.txt", std::ios::binary);
                if (!dst) {
                    std::cout << "Failed to open destination file" << std::endl;
                    return false;
                }
                dst << src.rdbuf();
            }
            
            form->SimulateCommand();  // Call the public method
            std::cout << "Processed test file: " << testFiles[i] << std::endl;
        }
    }
    catch (System::Exception^ ex) {
        System::Console::WriteLine("Test failed with error: {0}", ex->Message);
        success = false;
    }
    catch (std::exception& ex) {
        std::cout << "Test failed with error: " << ex.what() << std::endl;
        success = false;
    }
    catch (...) {
        std::cout << "Test failed with unknown error" << std::endl;
        success = false;
    }

    return success;
}