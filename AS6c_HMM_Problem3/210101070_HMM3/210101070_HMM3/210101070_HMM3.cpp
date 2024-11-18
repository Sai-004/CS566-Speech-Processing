// 210101070_HMM3.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iomanip>

using namespace std;

// Load matrix from file
void loadMatrixFromFile(const string& filepath, double** mat, int numRows, int numCols) 
{
    ifstream inputFile(filepath.c_str());
    if (!inputFile) {
        cerr << "Error: Unable to open the file " << filepath.c_str() << endl;
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            inputFile >> mat[i][j];
        }
    }
    inputFile.close();
}

// Load array from file
void loadArrayFromFile(const string& filepath, double* vec, int length) 
{
    ifstream inputFile(filepath.c_str());
    if (!inputFile) {
        cerr << "Error: Unable to open the file " << filepath.c_str() << endl;
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < length; ++i) {
        inputFile >> vec[i];
    }
    inputFile.close();
}

// For printing matrix for debugging
void displayMatrix(double** mat, int rows, int cols, const string& matrixName) 
{
    cout << "Displaying Matrix: " << matrixName.c_str() << endl;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            cout << fixed << setprecision(4) << mat[i][j] << "\t";
        }
        cout << endl;
    }
    cout << endl;
}

// Calculate sequence probability using forward variables (Alpha)
double calculateSequenceProbability(double** Alpha, int T, int numStates) {
    double probability = 0.0;
    for (int i = 0; i < numStates; ++i) {
        probability += Alpha[T-1][i];
    }
    return probability;
}

// Baum-Welch algorithm to estimate the HMM parameters
void baumWelch(double** A, double** B, double* Pi, int* O, int T, int numStates, int numObservations, int maxIterations = 100) 
{
    // Allocate memory for temporary matrices
    double** Alpha = new double*[T];
    double** Beta = new double*[T];
    double** Gamma = new double*[T];
    double*** Xi = new double**[T-1];
    
    for (int t = 0; t < T; t++) {
        Alpha[t] = new double[numStates]();
        Beta[t] = new double[numStates]();
        Gamma[t] = new double[numStates]();
        if (t < T-1) {
            Xi[t] = new double*[numStates];
            for (int i = 0; i < numStates; i++) {
                Xi[t][i] = new double[numStates]();
            }
        }
    }

    double previousProbability = 0.0;

    // Start iteration loop for Baum-Welch
    for (int iter = 0; iter < maxIterations; ++iter) {
        // Step 1: Forward Algorithm (Alpha)
        // Initialization
        for (int i = 0; i < numStates; ++i) {
            Alpha[0][i] = Pi[i] * B[i][O[0]];
        }
        // Induction
        for (int t = 1; t < T; ++t) {
            for (int j = 0; j < numStates; ++j) {
                Alpha[t][j] = 0.0;
                for (int i = 0; i < numStates; ++i) {
                    Alpha[t][j] += Alpha[t - 1][i] * A[i][j];
                }
                Alpha[t][j] *= B[j][O[t]];
            }
        }

        // Calculate probability for this iteration
        double currentProbability = calculateSequenceProbability(Alpha, T, numStates);
        double probImprovement = fabs(currentProbability - previousProbability);

        cout << "Iteration: " << iter + 1 << endl;
        cout << "Sequence Probability: " << scientific << setprecision(10) << currentProbability << endl;
        if (iter > 0) {
            cout << "Probability Improvement: " << scientific << setprecision(10) << probImprovement << endl;
        }
        cout << endl;

        previousProbability = currentProbability;

        // Step 2: Backward Algorithm (Beta)
        // Initialization
        for (int i = 0; i < numStates; ++i) {
            Beta[T - 1][i] = 1.0;
        }
        // Induction
        for (int t = T - 2; t >= 0; --t) {
            for (int i = 0; i < numStates; ++i) {
                Beta[t][i] = 0.0;
                for (int j = 0; j < numStates; ++j) {
                    Beta[t][i] += A[i][j] * B[j][O[t + 1]] * Beta[t + 1][j];
                }
            }
        }

        // Step 3: Compute Gamma and Xi
        for (int t = 0; t < T - 1; ++t) {
            double sum = 0.0;
            for (int i = 0; i < numStates; ++i) {
                for (int j = 0; j < numStates; ++j) {
                    Xi[t][i][j] = Alpha[t][i] * A[i][j] * B[j][O[t + 1]] * Beta[t + 1][j];
                    sum += Xi[t][i][j];
                }
            }
            // Normalize Xi and calculate Gamma
            for (int i = 0; i < numStates; ++i) {
                Gamma[t][i] = 0.0;
                for (int j = 0; j < numStates; ++j) {
                    Xi[t][i][j] /= sum;
                    Gamma[t][i] += Xi[t][i][j];
                }
            }
        }

        // Special case for Gamma_T-1
        double sumGamma = 0.0;
        for (int i = 0; i < numStates; ++i) {
            Gamma[T - 1][i] = Alpha[T - 1][i];
            sumGamma += Gamma[T - 1][i];
        }
        for (int i = 0; i < numStates; ++i) {
            Gamma[T - 1][i] /= sumGamma;
        }

        // Step 4: Re-estimate A, B, and Pi
        // Re-estimate Pi
        for (int i = 0; i < numStates; ++i) {
            Pi[i] = Gamma[0][i];
        }

        // Re-estimate A
        for (int i = 0; i < numStates; ++i) {
            for (int j = 0; j < numStates; ++j) {
                double numerator = 0.0;
                double denominator = 0.0;
                for (int t = 0; t < T - 1; ++t) {
                    numerator += Xi[t][i][j];
                    denominator += Gamma[t][i];
                }
                A[i][j] = numerator / denominator;
            }
        }

        // Re-estimate B
        for (int i = 0; i < numStates; ++i) {
            for (int k = 0; k < numObservations; ++k) {
                double numerator = 0.0;
                double denominator = 0.0;
                for (int t = 0; t < T; ++t) {
                    if (O[t] == k) {
                        numerator += Gamma[t][i];
                    }
                    denominator += Gamma[t][i];
                }
                B[i][k] = numerator / denominator;
            }
        }

        // Optional: Add convergence check
        if (iter > 0 && probImprovement == 0) {
            cout << "Converged at iteration " << iter + 1 << endl;
            break;
        }

        // Display updated matrices A, B, and Pi
        displayMatrix(A, numStates, numStates, "Updated A Matrix");
        displayMatrix(B, numStates, numObservations, "Updated B Matrix");
        cout << "Updated Pi Vector: ";
        for (int i = 0; i < numStates; ++i) {
            cout << fixed << setprecision(4) << Pi[i] << "\t";
        }
        cout << endl << endl;
    }

    // Clean up memory
    for (int t = 0; t < T; t++) {
        delete[] Alpha[t];
        delete[] Beta[t];
        delete[] Gamma[t];
        if (t < T-1) {
            for (int i = 0; i < numStates; i++) {
                delete[] Xi[t][i];
            }
            delete[] Xi[t];
        }
    }
    delete[] Alpha;
    delete[] Beta;
    delete[] Gamma;
    delete[] Xi;
}

int main() 
{
    const int numStates = 5;
    const int numObservations = 32;
    
    // Allocate memory for matrices and arrays
    double** A = new double*[numStates];
    double** B = new double*[numStates];
    double* Pi = new double[numStates];
    
    for (int i = 0; i < numStates; i++) {
        A[i] = new double[numStates];
        B[i] = new double[numObservations];
    }

    // Load matrices and arrays from files
    loadMatrixFromFile("files/A_ijMatrix.txt", A, numStates, numStates);
    loadMatrixFromFile("files/B_jkMatrix.txt", B, numStates, numObservations);
    loadArrayFromFile("files/P_iVector.txt", Pi, numStates);

    // Read observation sequence from file
    ifstream observationFile("files/ObservationSequence.txt");
    if (!observationFile) {
        cerr << "Error: Unable to open observation sequence file" << endl;
        return EXIT_FAILURE;
    }

    // First count the number of observations
    int T = 0;
    int temp;
    while (observationFile >> temp) {
        T++;
    }
    observationFile.clear();
    observationFile.seekg(0);

    // Allocate and read observation sequence
    int* O = new int[T];
    for (int i = 0; i < T; i++) {
        observationFile >> temp;
        O[i] = temp - 1;  // Convert to zero-based index
    }
    observationFile.close();

    // Call the Baum-Welch algorithm
    baumWelch(A, B, Pi, O, T, numStates, numObservations);

    // Clean up memory
    for (int i = 0; i < numStates; i++) {
        delete[] A[i];
        delete[] B[i];
    }
    delete[] A;
    delete[] B;
    delete[] Pi;
    delete[] O;

    return 0;
}