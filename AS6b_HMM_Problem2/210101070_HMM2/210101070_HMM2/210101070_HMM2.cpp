// 210101070_HMM2.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <algorithm>

using namespace std;

// Load matrix data from file
void loadMatrixFromFile(const string &filepath, vector<vector<double>> &mat, int numRows, int numCols)
{
	ifstream inputFile(filepath.c_str());
	if (!inputFile)
	{
		cerr << "Error: Unable to open the file " << filepath.c_str() << endl;
		exit(EXIT_FAILURE); // Exit on error
	}

	mat.resize(numRows, vector<double>(numCols));
	for (int i = 0; i < numRows; ++i)
	{
		for (int j = 0; j < numCols; ++j)
		{
			inputFile >> mat[i][j];
		}
	}
	inputFile.close();
}

// Load vector data from file
void loadVectorFromFile(const string &filepath, vector<double> &vec, int length)
{
	ifstream inputFile(filepath.c_str());
	if (!inputFile)
	{
		cerr << "Error: Unable to open the file " << filepath.c_str() << endl;
		exit(EXIT_FAILURE);
	}

	vec.resize(length);
	for (int i = 0; i < length; ++i)
	{
		inputFile >> vec[i];
	}
	inputFile.close();
}

// Viterbi algorithm to determine the most likely sequence of hidden states
vector<int> ViterbiAlgorithm(const vector<vector<double>> &A, const vector<vector<double>> &B, const vector<double> &init_Pi, const vector<int> &ObservationSequence)
{
	int numStates = A.size();						 // Number of states (N)
	int sequenceLength = ObservationSequence.size(); // Observation sequence length (T)

	// Table to store the probability of each state at each time step
	vector<vector<double>> DeltaTable(sequenceLength, vector<double>(numStates, 0.0));

	// Path table to track back-pointers for state transitions
	vector<vector<int>> pathTable(sequenceLength, vector<int>(numStates, 0));

	// Initializing process
	for (int i = 0; i < numStates; i++)
	{
		DeltaTable[0][i] = init_Pi[i] * B[i][ObservationSequence[0]];
		pathTable[0][i] = 0;
	}

	// Recursion process
	for (int t = 1; t < sequenceLength; t++)
	{
		for (int j = 0; j < numStates; j++)
		{
			double maxProb = -1.0;
			int bestState = 0;
			for (int i = 0; i < numStates; i++)
			{
				double prob = DeltaTable[t - 1][i] * A[i][j];
				if (prob > maxProb)
				{
					maxProb = prob;
					bestState = i;
				}
			}
			DeltaTable[t][j] = maxProb * B[j][ObservationSequence[t]];
			pathTable[t][j] = bestState;
		}
	}

	// Termination process: Find the state with the highest probability at the last time step
	double maxProb = -1.0;
	int finalState = 0;
	for (int i = 0; i < numStates; i++)
	{
		if (DeltaTable[sequenceLength - 1][i] > maxProb)
		{
			maxProb = DeltaTable[sequenceLength - 1][i];
			finalState = i;
		}
	}

	// Backtracking to reconstruct the most likely state sequence
	vector<int> stateSequence(sequenceLength);
	stateSequence[sequenceLength - 1] = finalState;

	for (int t = sequenceLength - 2; t >= 0; t--)
	{
		stateSequence[t] = pathTable[t + 1][stateSequence[t + 1]];
	}

	return stateSequence;
}

int main()
{
	vector<vector<double>> transitionMatrix;  // A: State transition probabilities
	vector<vector<double>> observationMatrix; // B: Observation likelihoods
	vector<double> initialProbabilities;	  // Pi: Initial state distribution
	vector<int> observationSequence;		  // O: Observation sequence

	// Load transition matrix, observation matrix, and initial probabilities from files
	loadMatrixFromFile("files/A_ijMatrix.txt", transitionMatrix, 5, 5);
	loadMatrixFromFile("files/B_jkMatrix.txt", observationMatrix, 5, 32);
	loadVectorFromFile("files/P_iVector.txt", initialProbabilities, 5);

	// Read observation sequence from file
	ifstream observationFile("files/ObservationSequence.txt");
	if (!observationFile)
	{
		cerr << "Error: Unable to open observation sequence file" << endl;
		return EXIT_FAILURE;
	}

	int obsValue;
	while (observationFile >> obsValue)
	{
		observationSequence.push_back(obsValue - 1); // Adjusting for zero-indexed sequence
	}
	observationFile.close();

	// Execute the Viterbi algorithm
	vector<int> mostLikelyStateSeq = ViterbiAlgorithm(transitionMatrix, observationMatrix, initialProbabilities, observationSequence);

	// Output the most probable sequence of hidden states
	cout << "Most probable state sequence: ";
	for (int i = 0; i < mostLikelyStateSeq.size(); ++i)
	{
		cout << mostLikelyStateSeq[i] + 1 << " "; // Adjusting back to one-indexed for output
	}
	cout << endl;

	return 0;
}