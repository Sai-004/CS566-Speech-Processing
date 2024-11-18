#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <filesystem>

// Custom weight configuration with a different distribution
const std::vector<double> customWeights = {
    1.2, 2.7, 6.5, 12.3, 18.1, 21.9,
    24.6, 32.1, 41.5, 49.8, 55.6, 60.2};

// Euclidean distance with weighted squared difference calculation
double calculateWeightedDistance(
    const std::vector<double> &vector1,
    const std::vector<double> &vector2,
    const std::vector<double> &weights)
{
    double totalDistance = 0.0;
    for (size_t i = 0; i < vector1.size(); ++i)
    {
        double diff = vector1[i] - vector2[i];
        totalDistance += weights[i] * diff * diff;
    }
    return std::sqrt(totalDistance);
}

// Alternative codebook loading with error handling
class CodebookManager
{
private:
    std::vector<std::vector<double>> codebookEntries;

public:
    bool loadCodebook(const std::string &filename)
    {
        std::ifstream file(filename);
        if (!file.is_open())
        {
            std::cerr << "Error: Unable to open codebook file " << filename << std::endl;
            return false;
        }

        codebookEntries.clear();
        std::string line;
        while (std::getline(file, line))
        {
            std::istringstream ss(line);
            std::vector<double> codeVector;
            double coefficient;

            while (ss >> coefficient)
            {
                codeVector.push_back(coefficient);
                if (ss.peek() == ',')
                    ss.ignore();
            }

            if (codeVector.size() == customWeights.size())
            {
                codebookEntries.push_back(codeVector);
            }
        }

        return !codebookEntries.empty();
    }

    int findNearestCodebookIndex(const std::vector<double> &inputVector) const
    {
        int nearestIndex = -1;
        double minDistance = std::numeric_limits<double>::max();

        for (size_t i = 0; i < codebookEntries.size(); ++i)
        {
            double distance = calculateWeightedDistance(
                inputVector,
                codebookEntries[i],
                customWeights);

            if (distance < minDistance)
            {
                minDistance = distance;
                nearestIndex = static_cast<int>(i);
            }
        }

        return nearestIndex;
    }

    size_t getCodebookSize() const
    {
        return codebookEntries.size();
    }
};

// Observation sequence generator with improved file handling
class ObservationGenerator
{
public:
    static bool generateObservationSequences(
        const std::string &inputFilename,
        const CodebookManager &codebookManager,
        const std::string &outputDir = "observation")
    {
        // Create output directory if it doesn't exist
        std::filesystem::create_directories(outputDir);

        std::ifstream inputFile(inputFilename);
        if (!inputFile.is_open())
        {
            std::cerr << "Error: Unable to open input file " << inputFilename << std::endl;
            return false;
        }

        std::string line, currentFilename;
        std::ofstream outputFile;

        while (std::getline(inputFile, line))
        {
            std::istringstream ss(line);
            std::string filename;
            int index;
            std::vector<double> cepstralCoeffs(customWeights.size());

            // Parse input line
            std::getline(ss, filename, ',');
            ss >> index;
            ss.ignore();

            for (size_t i = 0; i < cepstralCoeffs.size(); ++i)
            {
                ss >> cepstralCoeffs[i];
                if (ss.peek() == ',')
                    ss.ignore();
            }

            // Manage output file
            if (filename != currentFilename)
            {
                if (outputFile.is_open())
                    outputFile.close();
                currentFilename = filename;
                outputFile.open(outputDir + "/" + filename);
            }

            // Find and write observation
            int observationIndex = codebookManager.findNearestCodebookIndex(cepstralCoeffs);
            outputFile << observationIndex << ' ';
        }

        return true;
    }
};

int main()
{
    CodebookManager codebookManager;

    if (!codebookManager.loadCodebook("generatedCodeBook.csv"))
    {
        std::cerr << "Codebook loading failed." << std::endl;
        return 1;
    }

    if (!ObservationGenerator::generateObservationSequences(
            "CepstralCoefficients.csv",
            codebookManager))
    {
        std::cerr << "Observation sequence generation failed." << std::endl;
        return 1;
    }

    std::cout << "Observation sequences successfully generated." << std::endl;
    return 0;
}