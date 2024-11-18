#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <limits>
#include <string>
#include <algorithm>

// Constants
constexpr int DIMENSIONS = 12;
constexpr int CODEBOOK_SIZE = 32;
constexpr double EPSILON = 0.03;
constexpr double CONVERGENCE_THRESHOLD = 0.0001;
constexpr int MAX_ITERATIONS = 100;

// Global variables
std::vector<std::vector<double>> universe;
const std::vector<double> weights = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};

// Calculate Tokhura distance between two vectors
double calculateDistance(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    double distance = 0.0;
    for (int i = 0; i < DIMENSIONS; i++) {
        double diff = vec1[i] - vec2[i];
        distance += weights[i] * diff * diff;
    }
    return distance;
}

// Split codebook vectors
std::vector<std::vector<double>> splitCodebook(const std::vector<std::vector<double>>& codebook) {
    std::vector<std::vector<double>> newCodebook;
    
    for (const auto& codeVector : codebook) {
        std::vector<double> vectorPlus(DIMENSIONS);
        std::vector<double> vectorMinus(DIMENSIONS);

        for (int i = 0; i < DIMENSIONS; i++) {
            vectorPlus[i] = codeVector[i] * (1.0 + EPSILON);
            vectorMinus[i] = codeVector[i] * (1.0 - EPSILON);
        }

        newCodebook.push_back(vectorPlus);
        newCodebook.push_back(vectorMinus);
    }

    return newCodebook;
}

// Apply K-means clustering
double kMeans(std::vector<std::vector<double>>& codebook) {
    std::vector<int> clusterAssignments(universe.size());
    std::vector<int> clusterSizes(codebook.size());
    double prevDistortion = std::numeric_limits<double>::max();
    double distortion;
    int iteration = 0;

    do {
        distortion = 0.0;
        std::fill(clusterSizes.begin(), clusterSizes.end(), 0);
        std::vector<std::vector<double>> newCodebook(codebook.size(), std::vector<double>(DIMENSIONS, 0.0));

        // Assign vectors to nearest codebook entry
        for (size_t i = 0; i < universe.size(); i++) {
            double minDist = std::numeric_limits<double>::max();
            int bestMatch = 0;

            for (size_t j = 0; j < codebook.size(); j++) {
                double dist = calculateDistance(universe[i], codebook[j]);
                if (dist < minDist) {
                    minDist = dist;
                    bestMatch = j;
                }
            }

            clusterAssignments[i] = bestMatch;
            clusterSizes[bestMatch]++;
            distortion += minDist;

            for (int d = 0; d < DIMENSIONS; d++) {
                newCodebook[bestMatch][d] += universe[i][d];
            }
        }

        // Calculate new centroids
        for (size_t i = 0; i < codebook.size(); i++) {
            if (clusterSizes[i] > 0) {
                for (int d = 0; d < DIMENSIONS; d++) {
                    newCodebook[i][d] /= clusterSizes[i];
                }
            }
        }

        codebook = newCodebook;
        iteration++;

        double avgDistortion = distortion / universe.size();
        std::cout << "Iteration " << iteration << ": Average Distortion = " << avgDistortion << std::endl;

    } while (iteration < MAX_ITERATIONS && 
            std::abs(distortion - prevDistortion) > CONVERGENCE_THRESHOLD * prevDistortion);

    return distortion / universe.size();
}

// Read cepstral coefficients from CSV file
bool loadCepstralCoefficients(const std::string& filename) {
    std::ifstream file(filename.c_str());
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return false;
    }

    std::string line;
    // Skip header line
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string cell;
        std::vector<double> coefficients;

        // Skip filename and frame index
        std::getline(ss, cell, ',');  // Skip filename
        std::getline(ss, cell, ',');  // Skip frame index

        for (int i = 0; i < DIMENSIONS; i++) {
            std::getline(ss, cell, ',');
            coefficients.push_back(std::stod(cell));
        }

        universe.push_back(coefficients);
    }

    std::cout << "Loaded " << universe.size() << " vectors from the file." << std::endl;
    return true;
}

// Generate codebook using LBG algorithm
std::vector<std::vector<double>> generateCodebook() {
    if (universe.empty()) {
        std::cerr << "Error: No data loaded" << std::endl;
        return std::vector<std::vector<double>>();
    }

    // Start with centroid of entire universe
    std::vector<std::vector<double>> codebook(1, std::vector<double>(DIMENSIONS, 0.0));
    for (const auto& vector : universe) {
        for (int i = 0; i < DIMENSIONS; i++) {
            codebook[0][i] += vector[i];
        }
    }
    for (int i = 0; i < DIMENSIONS; i++) {
        codebook[0][i] /= universe.size();
    }

    // Split codebook until desired size is reached
    while (codebook.size() < CODEBOOK_SIZE) {
        std::cout << "\nSplitting codebook from " << codebook.size() 
                  << " to " << codebook.size() * 2 << " vectors..." << std::endl;
        
        codebook = splitCodebook(codebook);
        double distortion = kMeans(codebook);
        
        std::cout << "Final distortion for " << codebook.size() 
                  << " codebook vectors: " << distortion << std::endl;
    }

    return codebook;
}

// Save codebook to file
void saveCodebook(const std::vector<std::vector<double>>& codebook, const std::string& filename) {
    std::ofstream file(filename.c_str());
    if (!file.is_open()) {
        std::cerr << "Error: Could not open output file" << std::endl;
        return;
    }

    for (const auto& vector : codebook) {
        for (int i = 0; i < DIMENSIONS; i++) {
            file << vector[i];
            if (i < DIMENSIONS - 1) {
                file << ",";
            }
        }
        file << "\n";
    }

    std::cout << "Codebook saved to " << filename << std::endl;
}

int main() {
    // Load the cepstral coefficients
    if (!loadCepstralCoefficients("CepstralCoefficients.csv"))
    {
        return 1;
    }

    // Generate the codebook
    std::cout << "\nGenerating codebook...\n" << std::endl;
    std::vector<std::vector<double>> codebook = generateCodebook();
    std::cout<<universe.size();
    // Save the codebook
    saveCodebook(codebook, "generatedCodeBook.csv");

    return 0;
}