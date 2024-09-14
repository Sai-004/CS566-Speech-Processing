// Yes_No_Detection.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Constants for frame size and total data
#define FRAME_SIZE 400  // Adjusted for 25ms of phoneme identification duration (16kHz sampling rate -> 0.025s * 16000 = 400 samples/frame)
#define MAX_DATA_SIZE 16000  // Assumed maximum size of the audio data

// Calculate energy of a frame
double calculateEnergy(short* frame, int frameSize) {
    double energy = 0.0;
    for (int i = 0; i < frameSize; i++) {
        energy += frame[i] * frame[i];
    }
    return energy;
}

// Calculate Zero Crossing Rate (ZCR) of a frame
double calculateZCR(short* frame, int frameSize) {
    int zcr = 0;
    for (int i = 1; i < frameSize; i++) {
        if ((frame[i] > 0 && frame[i - 1] < 0) || (frame[i] < 0 && frame[i - 1] > 0)) {
            zcr++;
        }
    }
    return (double)zcr / frameSize;
}

// Load audio data from a text file into a static array
int loadAudioData(const char* fileName, short* audioData) {
    FILE* file = fopen(fileName, "r");
    if (file == NULL) {
        printf("Error opening file\n");
        return -1;
    }

    int size = 0;
    short sample;

    // Skip first 5 lines (coz of metadata, can be commented if no meta data)
    char buffer[100];
    for (int i = 0; i < 5; i++) {
        fgets(buffer, sizeof(buffer), file);
    }

    while (fscanf(file, "%hd", &sample) != EOF && size < MAX_DATA_SIZE) {
        audioData[size++] = sample;
    }

    fclose(file);
    return size;  // Return the number of samples loaded
}

// Process files for training
double train() {
    const char* filesYes[] = {
        "audio/yes_1.txt", "audio/yes_2.txt", "audio/yes_3.txt", "audio/yes_4.txt", "audio/yes_5.txt"
    };
    const char* filesNo[] = {
        "audio/no_1.txt", "audio/no_2.txt", "audio/no_3.txt", "audio/no_4.txt", "audio/no_5.txt"
    };

    double totalZCRYes = 0.0, totalEnergyYes = 0.0;
    double totalZCRNo = 0.0, totalEnergyNo = 0.0;
    int countYes = 0, countNo = 0;

    // Process 'Yes' files
    for (int i = 0; i < 5; i++) {
        short audioData[MAX_DATA_SIZE];
        int dataSize = loadAudioData(filesYes[i], audioData);
        if (dataSize == -1) continue;

        int frameCount = 0;
        double fileZCR = 0.0, fileEnergy = 0.0;

        for (int j = 0; j + FRAME_SIZE <= dataSize; j += FRAME_SIZE) {
            short frame[FRAME_SIZE];
            memcpy(frame, &audioData[j], FRAME_SIZE * sizeof(short));

            double zcr = calculateZCR(frame, FRAME_SIZE);
            double energy = calculateEnergy(frame, FRAME_SIZE);

            fileZCR += zcr;
            fileEnergy += energy;
            frameCount++;
        }

        totalZCRYes += fileZCR / frameCount;
        totalEnergyYes += fileEnergy / frameCount;
        countYes++;
    }

    // Process 'No' files
    for (int i = 0; i < 5; i++) {
        short audioData[MAX_DATA_SIZE];
        int dataSize = loadAudioData(filesNo[i], audioData);
        if (dataSize == -1) continue;

        int frameCount = 0;
        double fileZCR = 0.0, fileEnergy = 0.0;

        for (int j = 0; j + FRAME_SIZE <= dataSize; j += FRAME_SIZE) {
            short frame[FRAME_SIZE];
            memcpy(frame, &audioData[j], FRAME_SIZE * sizeof(short));

            double zcr = calculateZCR(frame, FRAME_SIZE);
            double energy = calculateEnergy(frame, FRAME_SIZE);

            fileZCR += zcr;
            fileEnergy += energy;
            frameCount++;
        }

        totalZCRNo += fileZCR / frameCount;
        totalEnergyNo += fileEnergy / frameCount;
        countNo++;
    }

    // Calculate average ZCR and Energy for 'Yes' and 'No' files
    double avgZCRYes = totalZCRYes / countYes;
    double avgEnergyYes = totalEnergyYes / countYes;
    double avgZCRNo = totalZCRNo / countNo;
    double avgEnergyNo = totalEnergyNo / countNo;

    // Calculate threshold ZCR (midpoint between 'Yes' and 'No' average ZCRs)
    double zcrThreshold = (avgZCRYes + avgZCRNo) / 2.0;
    double energyThreshold = (avgEnergyYes + avgEnergyNo) / 2.0;

    // Write the results to a CSV file
    FILE* csvFile = fopen("results/training_results.csv", "w");
    if (csvFile == NULL) {
        printf("Error creating CSV file\n");
        return 0.0;
    }
    fprintf(csvFile, "File,Average ZCR,Average Energy\n");
    fprintf(csvFile, "Yes Files,%.2f,%.2f\n", avgZCRYes, avgEnergyYes);
    fprintf(csvFile, "No Files,%.2f,%.2f\n", avgZCRNo, avgEnergyNo);
    fprintf(csvFile, "ZCR Threshold,%.2f\n", zcrThreshold);
    fprintf(csvFile, "Energy Threshold,%.2f\n", energyThreshold);
    fclose(csvFile);

    printf("Training complete. ZCR Threshold: %.2f, Energy Threshold: %.2f\n", zcrThreshold, energyThreshold);
    printf("Number of files processed: 10\n");

	return zcrThreshold;
}

// Process a file for testing
void test(double zcrThreshold) {
    const char* testFile = "test/no_test.txt"; // test files. Make sure to change to the required file name to test before running.

    short audioData[MAX_DATA_SIZE];
    int dataSize = loadAudioData(testFile, audioData);
    if (dataSize == -1) return;

    int frameCount = 0;
    double totalZCR = 0.0;

    for (int i = 0; i + FRAME_SIZE <= dataSize; i += FRAME_SIZE) {
        short frame[FRAME_SIZE];
        memcpy(frame, &audioData[i], FRAME_SIZE * sizeof(short));

        double zcr = calculateZCR(frame, FRAME_SIZE);
        totalZCR += zcr;
        frameCount++;
    }

    double avgZCR = totalZCR / frameCount;

    // Make decision based on average ZCR and threshold
    const char* decision = (avgZCR > zcrThreshold) ? "Yes" : "No";

    printf("Testing complete. Decision: %s\n", decision);
    printf("Number of files created: 1 (test_output.csv)\n");

    // Write the results to a CSV file
    FILE* csvFile = fopen("results/test_output.csv", "a");
    if (csvFile == NULL) {
        printf("Error creating CSV file\n");
        return;
    }
    fprintf(csvFile, "File,Average ZCR,Decision\n");
    fprintf(csvFile, "%s,%.2f,%s\n", testFile, avgZCR, decision);
    fclose(csvFile);
}

int _tmain(int argc, _TCHAR* argv[]) {
    int choice;
    double zcrThreshold = 0.0;  // Initial ZCR threshold

    while (1) {
        printf("Select an option:\n");
        printf("1. Process sample audio files and get ZCR Threshold\n");
        printf("2. Test an audio and Decide 'Yes' or 'No'\n");
        printf("3. Exit\n");
        printf("Enter your choice: ");
        scanf("%d", &choice);

        switch (choice) {
        case 1:
            zcrThreshold = train();  // update ZCR threshold
            break;
        case 2:
            if (zcrThreshold == 0.0) {
                printf("Please process the sample audio files first to determine the ZCR threshold.\n");
            } else {
                test(zcrThreshold);
            }
            break;
        case 3:
            printf("Exiting...\n");
            return 0;
        default:
            printf("Invalid choice. Please try again.\n");
        }
        printf("Done.\n\n");
    }

    return 0;
}