// KMeansAlgotrithm_Algorithm.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <cstdlib>
using namespace std;

// Define constants
#define DELTA 0.0001    // Convergence threshold
#define P 12            // Order of cepstral coefficients

const double tokhuraWt[P] = {1, 3, 7, 13, 19, 22, 25, 33, 42, 50, 56, 61};  // Tokhura weights

double ComputeTokhuraDistance(const vector<double> &p1, const vector<double> &p2)
{
    double dist = 0.0;
    for (int k = 0; k < P; ++k)
    {
		double diff = (p1[k]-p2[k]);
        dist += tokhuraWt[k] * diff * diff;
    }
    return dist;
}

int main()
{
    vector<vector<double>> unv;
    string line;
    ifstream file("universe.csv");
    
    // Reading universe data from CSV file
    while (getline(file, line))
    {
        vector<double> row(P);
        stringstream ss(line);
        for (int i = 0; i < P; ++i)
        {
            ss >> row[i];
            if (ss.peek() == ',')
                ss.ignore();
        }
        unv.push_back(row);
    }

    srand(time(0)); // Random seed

    // K-Means algorithm
    int M = unv.size(); // Size of the universe
    int k = 8;          // Number of clusters

    vector<vector<double>> CB(k, vector<double>(P)); // Codebook
    vector<int> assign(M, 0);                        // Region assignment for each vector
    vector<double> distortion;                       // Distortion values in each iteration

    // Initializing codebook with random universe vectors
    for (int i = 0; i < k; ++i)
        CB[i] = unv[rand() % M];

    int itr = 0;
    double curr_distortion = 0.0;
    double prev_distortion = numeric_limits<double>::max();

    while (abs(prev_distortion - curr_distortion) > DELTA)
    {
        prev_distortion = curr_distortion;
        curr_distortion = 0.0;

        // Applying Nearest neighbor rule: Assign each vector to nearest cluster
        for (int i = 0; i < M; ++i)
        {
            double min_dist = numeric_limits<double>::max();
            int nearest_cluster = 0;
            for (int j = 0; j < k; ++j)
            {
                double dist = ComputeTokhuraDistance(unv[i], CB[j]);
                if (dist < min_dist)
                {
                    min_dist = dist;
                    nearest_cluster = j;
                }
            }
            assign[i] = nearest_cluster;
            curr_distortion += min_dist;
        }

        // Calculating average distortion
        double avg_distortion = curr_distortion / M;
        cout << "Iteration " << itr << "\t--> Average Distortion = " << avg_distortion << endl;

        // Updating centroids for each region
        vector<vector<double>> new_CB(k, vector<double>(P, 0.0));
        vector<int> cluster_sizes(k, 0);

        for (int i = 0; i < M; ++i)
        {
            int cluster = assign[i];
            for (int j = 0; j < P; ++j)
            {
                new_CB[cluster][j] += unv[i][j];
            }
            cluster_sizes[cluster]++;
        }

        for (int i = 0; i < k; ++i)
        {
            if (cluster_sizes[i] > 0)
            {
                for (int j = 0; j < P; ++j)
                {
                    new_CB[i][j] /= cluster_sizes[i];
                }
            }
            else
                new_CB[i] = unv[rand() % M]; // assigning random vectors to empty clusters
        }

        CB = new_CB;
        itr++;
    }

    // Printing codebook vectors
    cout << "\nFinal codebook vectors:" << endl;
	printf("-----------------------------------------------------------------------------------------------------------------------------------------------\n");
    for (int i = 0; i < k; ++i)
    {
        for (int j = 0; j < P; ++j)
        {
			printf("| %9.6lf ", CB[i][j]);
        }
        printf("|\n");
    }
	printf("------------------------------------------------------------------------------------------------------------------------------------------------\n");

    system("pause");
    return 0;
}