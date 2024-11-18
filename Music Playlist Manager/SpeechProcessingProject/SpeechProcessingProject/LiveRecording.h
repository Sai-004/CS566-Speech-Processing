#ifndef LIVE_RECORDING_H
#define LIVE_RECORDING_H

#include "stdafx.h"
#include <vector>
#include <string>
#include <Windows.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define FRAME_SIZE 320
#define HOP_LENGTH 240
#define P 12
#define SAMPLE_RATE 16025
#define ARRAY_SIZE (int)(SAMPLE_RATE * 1.2 + 1)

struct Frame {
    std::vector<double> data;
    std::vector<double> R;
    std::vector<double> A;
    std::vector<double> C;
};

class LiveRecording {
private:
    short int waveIn[ARRAY_SIZE];
    std::vector<double> weightsArray;
    Frame frame;

    void ListInputDevices();
    std::vector<double> readLiveData();
    void normalizeData(std::vector<double>& data);
    void applyHammingWindow(std::vector<double>& frame);
    std::vector<double> computeAutocorrelation(const std::vector<double>& frame, int p);
    std::vector<double> levinsonDurbin(const std::vector<double>& R);
    std::vector<double> computeCepstral(const std::vector<double>& A);
    void applyRaisedSineWindow(std::vector<double>& C);
    double tokhuraDistance(const std::vector<double>& c1, const std::vector<double>& c2);
    int findClosestCodebookIndex(const std::vector<double>& cepstralCoeffs, 
                                const std::vector<std::vector<double>>& codebook);
    std::vector<std::vector<double>> loadCodebook(const std::string& filename);

public:
    LiveRecording();  // Constructor declaration
    void StartRecording();
    void SaveToFile();
    void ProcessRecording();
    void GenerateObservationSequence();
    std::string GetObservationFilePath() const { return "live_observation.txt"; }
};

#endif