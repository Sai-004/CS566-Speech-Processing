// VowelRecognition : Entry point for the console application

#include "stdafx.h"
#include <stdio.h>
#include <cstring>
#include <limits>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

// Constants
#define FRAME_SIZE 320      // total sample in one frame
#define P 12                // Order of LPC
#define STABLE_FRAMES 5     // count of stable frames
#define PI 3.142857142857   // Value of pi for calculation
#define TRAINDATA_SIZE 20	// No.of train files/vowel
const double amplitudeThreshold = 5000.0; // Threshold multiplier for amplitude normalization

/* Global Variables */
// For max sample value, dc shift and normalization factor
double maxAmplitude, dcShift, nFactor;

// To store Samples, Energy, Frames, Tokhura distance
double Samples[100000], energyValues[100000], steadyFrames[STABLE_FRAMES][FRAME_SIZE], tokhuraDistance[5];

// To store framwise Ri, Ai, Ci, avgerage Ci, temp Ci, all Ci, reference Ci
double R[STABLE_FRAMES][P + 1], A[STABLE_FRAMES][P + 1], C[STABLE_FRAMES][P + 1];
double avgCi[25][P + 1], Ci[STABLE_FRAMES][P + 1], allCi[100][STABLE_FRAMES][P + 1], restoredCi[STABLE_FRAMES][P + 1];

// Tokhura distance calculation weights
double tokhuraWeights[] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};

// Indexing helpers to store values of ci to allCi
int indexCoeff = 0, fileCount = 0;

// Vowel set
char vowels[5] = {'a', 'e', 'i', 'o', 'u'};

// Sample array size, Energy array size, start and end points of steady frames
long int sampleSize, energySize, start_i, end_i;

// To find accuracy
int totalCorrect = 0, individualCorrect = 0;

// Apply raised Sine Window to Ci for each frame
void applySineWindow()
{
    long double sum = 0;
    for (int f = 0; f < STABLE_FRAMES; f++)
    {
        for (int coeff = 1; coeff <= P; coeff++)
        {
            sum = (P / 2) * sin((PI * coeff) / P);
            C[f][coeff] *= sum;
        }
    }
}

// Applying Hamming Window for each stable frame
void applyHammingWindow()
{
    for (int i = 0; i < STABLE_FRAMES; ++i)
    {
        for (int j = 0; j < FRAME_SIZE; ++j)
        {
            steadyFrames[i][j] *= 0.54 - 0.46 * cos(2 * PI * j / FRAME_SIZE - 1);
        }
    }
}

// Save Ci of each frame in allCi
void saveCis()
{
    // saving ci values to 3D matrix allCi
    for (int f = 0; f < STABLE_FRAMES; f++)
    {
        for (indexCoeff = 0; indexCoeff < P; indexCoeff++)
        {
            allCi[fileCount][f][indexCoeff + 1] = C[f][indexCoeff + 1];
        }
    }
    fileCount++;
}

// Store avgerage ci values in a file
void saveAvgCi()
{
    FILE *outputFile;
    char filename[80];
    int index = 0;

    // Save file of each vowel
    for (int vowel = 0; vowel < 5; vowel++)
    {
        sprintf(filename, "Reference/reference_ci_%c.txt", vowels[vowel]);
        outputFile = fopen(filename, "w");
        for (int f = 0; f < STABLE_FRAMES; f++)
        {
            for (int coeff = 0; coeff < P; coeff++)
            {
                // Calculate the average Ci for each file
                double sum = 0;
                for (int file = vowel * TRAINDATA_SIZE; file < (vowel + 1) * TRAINDATA_SIZE; file++)
                {
                    sum += allCi[file][f][coeff + 1];
                }
                sum /= TRAINDATA_SIZE;
                avgCi[index][coeff + 1] = sum;
                fprintf(outputFile, "%lf ", sum);
            }
            index++;
            fprintf(outputFile, "\n");
        }
        printf(">> Created %s\n", filename);
        fclose(outputFile);
    }
}

// Compute Cepstral coefficients Ci's
void ComputeCi()
{
    double sum = 0;

    for (int f = 0; f < STABLE_FRAMES; f++)
    {
        C[f][0] = log(R[f][0] * R[f][0]);

        for (int m = 1; m <= P; m++)
        {
            sum = 0;
            for (int k = 1; k < m; k++)
            {
                sum += (k * C[f][k] * A[f][m - k]) / (m * 1.0);
            }
            C[f][m] = A[f][m] + sum;
        }
    }

    applySineWindow(); // applying raised sine window on Ci's
    saveCis();         // save Ci's to allCi
}

// Find Ai's using Levinson Durbin Algorithm
void ComputeAi()
{
    double Alpha[P + 1][P + 1]; // 2D array for LPC coefficients
    double E[P + 1];            // Error predictions
    double K[P + 1];            // Reflection coefficients
    double sum = 0;

    for (int f = 0; f < STABLE_FRAMES; f++)
    {
        E[0] = R[f][0];
        for (int i = 1; i <= P; i++)
        {
            sum = 0;
            // Compute the prediction gain
            for (int j = 1; j <= i - 1; j++)
            {
                sum += Alpha[i - 1][j] * R[f][i - j];
            }

            K[i] = (R[f][i] - sum) / E[i - 1];
            Alpha[i][i] = K[i]; // alpha(i)_i = K_i

            for (int j = 1; j <= i - 1; j++)
            {
                Alpha[i][j] = Alpha[i - 1][j] - K[i] * Alpha[i - 1][i - j];
            }

            E[i] = (1 - (K[i] * K[i])) * E[i - 1];
        }

        // store LPC coefficients Ai's
        for (int i = 1; i <= P; i++)
        {
            A[f][i] = Alpha[P][i];
        }
    }

    // Compute Ci's
    ComputeCi();
}

// Computing autocorrelation coefficients Ri
void ComputeRi()
{
    for (int f = 0; f < STABLE_FRAMES; f++)
    {
        for (int m = 0; m <= P; m++)
        {
            R[f][m] = 0;
            for (int k = 0; k < FRAME_SIZE - m; k++)
            {
                R[f][m] += steadyFrames[f][k] * steadyFrames[f][k + m];
            }
        }
    }

    // Levinson Durbin Algorithm to calculate Ai's
    ComputeAi();
}

// Performing DC shift normalization
double getDCShift(char *filename)
{
    long int sample_count = 0;
    FILE *fp;
    char line[80];

    fp = fopen(filename, "r");

    if (fp == NULL)
    {
        printf("File not found\n");
        exit(1);
    }

    dcShift = 0;
    while (!feof(fp))
    {
        fgets(line, 80, fp);
        dcShift += atof(line);
        sample_count++;
    }
    dcShift /= sample_count; // Compute avg DC Shift

    fclose(fp);
    return dcShift;
}

// Initialize global variables from the input file
void setGlobalVariables(char *filename)
{
    FILE *fp;
    long int totalSample = 0;
    char line[100];

    fp = fopen(filename, "r");
    if (fp == NULL)
    {
        printf("Error opening file\n");
    }

	// Performing Amplitude normalization
    maxAmplitude = 0;
    while (!feof(fp))
    {
        fgets(line, 100, fp);
        if (!isalpha(line[0]))
        {
            totalSample++;
            if (maxAmplitude < abs(atoi(line)))
                maxAmplitude = abs(atoi(line));
        }
    }

	// Set normalization factor
    nFactor = (double)amplitudeThreshold / maxAmplitude;

	// Get DC Shift/Offset
    dcShift = getDCShift(filename);
    fclose(fp);
}

// To find stable frames
void findSteadyFrames()
{
    long int totalSample = 0, maxEnergyIndex = 0;
    int n = 0;
    double energy = 0, maxEnergy = 0;
    energySize = 0;

    // Find highest energy frame
	while (totalSample < sampleSize)
    {
        if (n == FRAME_SIZE)
        {
            // taking average
            energy /= FRAME_SIZE;

            if (maxEnergy < energy)
                maxEnergy = energy, maxEnergyIndex = energySize;

            energyValues[energySize++] = energy; // Store energy value

			// reset energy value and n for the next frame
            energy = 0;
            n = 0;
        }
        energy += pow(Samples[totalSample], 2);
        totalSample++;
        n++;
    }

    // Find start and end frames based on maximum energy
    start_i = maxEnergyIndex > 2 ? (maxEnergyIndex - 2) * FRAME_SIZE : 0;
    end_i = maxEnergyIndex < energySize - 3 ? (maxEnergyIndex + 3) * FRAME_SIZE : energySize * FRAME_SIZE;

    int f = 0;
    for (int i = start_i, j = 0; i < end_i; i++)
    {
        steadyFrames[f][j++] = Samples[i];
        if (j == FRAME_SIZE)
            f++, j = 0;
    }
}

// Process input file
void processFile(char *filename)
{
    char line[70];
    FILE *ip;

	// Applying normalizations and setting up the global variables
    setGlobalVariables(filename);

    ip = fopen(filename, "r");

    if (ip == NULL)
        printf("Error in opening file %s\n", filename);

    sampleSize = 0;
    while (!feof(ip))
    {
        fgets(line, 100, ip);

		// Skipping meta data, if any
        if (!isalpha(line[0]))
        {
            int y = atof(line);
            double normalizedSample = floor((y - dcShift) * nFactor);
            Samples[sampleSize++] = normalizedSample; // Store normalized sample
        }
    }

    fclose(ip);
    findSteadyFrames();   // Identify stable frames
    applyHammingWindow(); // Apply Hamming Window on stable frames
    ComputeRi();          // Compute Ri's on stable frames
}

// Train with sample data of each vowel
void train()
{
    printf("Training initiated...\n");
	char filename[100];

    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < TRAINDATA_SIZE; j++)
        {
            sprintf(filename, "TrainData/210101070_%c_%d.txt", vowels[i], j + 1);
			processFile(filename); // Process each training file
        }
		printf("Vowel '%c' > training done...\n",vowels[i]);
    }
    printf("Training completed.\n\n");
}

// Calculate the tokhura distance between Ci's
double calculateTokhuraDistance(FILE *ip)
{
    char line[3000];
    int f = 0;

	// read the reference file and store in restoredCi
    while (!feof(ip) && f < STABLE_FRAMES)
    {
        fgets(line, sizeof(line), ip);
        // Read cepstral coefficients from file
        sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &restoredCi[f][1], &restoredCi[f][2], &restoredCi[f][3],
               &restoredCi[f][4], &restoredCi[f][5], &restoredCi[f][6],
               &restoredCi[f][7], &restoredCi[f][8], &restoredCi[f][9],
               &restoredCi[f][10], &restoredCi[f][11], &restoredCi[f][12]);
        f++;
    }

    double totalDist = 0;
    for (int i = 0; i < STABLE_FRAMES; i++)
    {
        double frameDist = 0;
        for (int p = 1; p <= P; p++)
        {
            double d = (C[i][p] - restoredCi[i][p]);
            frameDist += tokhuraWeights[p - 1] * d * d; // Tokhura weighted distance
        }
        totalDist += frameDist / (P * 1.0); // Total distance
    }
    return totalDist / (STABLE_FRAMES * 1.0); // Average distance over all frames
}

// Predict vowel from test data using tokhura distance
char predictVowel()
{
    char filename[100];
    FILE *ip;
    double minDistance = DBL_MAX; // Initialize with maximum possible value
    char predictedVowel;

    for (int i = 0; i < 5; i++)
    {
        sprintf(filename, "Reference/reference_ci_%c.txt", vowels[i]);
        ip = fopen(filename, "r");
        if (ip == NULL)
            printf("Error in opening file: %s\n", filename);
		
		// Calculate distance from vowel's reference file
		double distance = calculateTokhuraDistance(ip);
        tokhuraDistance[i] = distance;
        
		if (minDistance > distance)
        {
            minDistance = distance;
            predictedVowel = vowels[i]; // Predict vowel by minimum distance
        }
    }
    return predictedVowel;
}

// Predict the vowel from test data
void test()
{
    char filename[100];
    fileCount = 0;

    int choice;
    printf("\n\nProcessing the test files...\n");

    for (int v = 0; v < 5; v++)
    {
        individualCorrect = 0;
        for (int file = TRAINDATA_SIZE; file < TRAINDATA_SIZE + 10; file++)
        {
            printf("\nTestData/210101070_%c_%d.txt", vowels[v], file + 1);
            sprintf(filename, "TestData/210101070_%c_%d.txt", vowels[v], file + 1);
            
			// Process each test file
			processFile(filename);

			// Predict the vowel
            char prediction = predictVowel();
            
			printf("\nTokhura Distance for {a, e, i, o, u} is {%lf, %lf, %lf, %lf, %lf}",
                   tokhuraDistance[0], tokhuraDistance[1], tokhuraDistance[2], tokhuraDistance[3], tokhuraDistance[4]);
			printf("\nFile: %s , vowel predicted ==> ", filename);
			printf("'%c'\n\n", prediction);

            if (prediction == vowels[v])
                totalCorrect++, individualCorrect++; // Track accuracy
        }
        printf("----------------- Accuracy for vowel %c is %.2lf %% -----------------\n\n", vowels[v], (individualCorrect / 10.0) * 100);
    }
    printf("---------------------------------------------------------------------------\n");
    printf("Overall accuracy is %.2lf %% \n", (totalCorrect / 50.0) * 100); // Print overall accuracy
    printf("---------------------------------------------------------------------------\n\n");
}

// Main function of the application
int main()
{
    train();	 // Train using 100 recordings of text data with 20 recordings per vowel
    saveAvgCi(); // Generate reference files
    test();      // Predict vowels in the test files

    printf("Testing completed.\n");
    system("pause");
    return 0;
}
