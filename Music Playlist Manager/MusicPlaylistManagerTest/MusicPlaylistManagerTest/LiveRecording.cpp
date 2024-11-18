#include "stdafx.h"
#include "LiveRecording.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
using namespace std;

#pragma comment(lib, "winmm.lib")
#define RECORD_SECONDS 1.2

LiveRecording::LiveRecording() {
    // Initialize weightsArray using push_back
    weightsArray.push_back(1.0);
    weightsArray.push_back(3.0);
    weightsArray.push_back(7.0);
    weightsArray.push_back(13.0);
    weightsArray.push_back(19.0);
    weightsArray.push_back(22.0);
    weightsArray.push_back(25.0);
    weightsArray.push_back(33.0);
    weightsArray.push_back(42.0);
    weightsArray.push_back(50.0);
    weightsArray.push_back(56.0);
    weightsArray.push_back(61.0);
}

void LiveRecording::ListInputDevices() {
    UINT numDevices = waveInGetNumDevs();
    WAVEINCAPS wic;
    for (UINT i = 0; i < numDevices; i++) {
        if (waveInGetDevCaps(i, &wic, sizeof(WAVEINCAPS)) == MMSYSERR_NOERROR) {
            std::cout << "Device ID " << i << ": " << wic.szPname << std::endl;
        }
    }
}

std::vector<std::vector<double>> LiveRecording::loadCodebook(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<std::vector<double>> codebook;
    std::string line;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::vector<double> codeVector;
        double coeff;
        while (ss >> coeff) {
            codeVector.push_back(coeff);
            if (ss.peek() == ',') ss.ignore();
        }
        codebook.push_back(codeVector);
    }
    return codebook;
}

void LiveRecording::StartRecording() {
    ListInputDevices();  // Show available input devices

    HWAVEIN hWaveIn;
    MMRESULT result;
    WAVEFORMATEX pFormat;

    pFormat.wFormatTag = WAVE_FORMAT_PCM;
    pFormat.nChannels = 1;                     // Mono audio
    pFormat.nSamplesPerSec = SAMPLE_RATE;      // Sample rate
    
	#ifdef BIT_DEPTH_16
		pFormat.wBitsPerSample = 16;               // 16-bit audio
		pFormat.nBlockAlign = sizeof(short int);   // Block alignment for 16-bit samples
		pFormat.nAvgBytesPerSec = SAMPLE_RATE * sizeof(short int);
	#else
		pFormat.wBitsPerSample = 8;                // 8-bit audio
		pFormat.nBlockAlign = 1;                   // Block alignment for 8-bit samples
		pFormat.nAvgBytesPerSec = SAMPLE_RATE;
	#endif

		UINT deviceID = 0;  // Set deviceID based on ListInputDevices output

    // Open the waveform audio input device
    result = waveInOpen(&hWaveIn, deviceID, &pFormat, 0L, 0L, CALLBACK_NULL);
    if (result != MMSYSERR_NOERROR) {
        cout << "Error opening waveform input device." << endl;
        return;
    }

    WAVEHDR WaveInHdr;
    WaveInHdr.lpData = (LPSTR)waveIn;  // Buffer for audio data
    WaveInHdr.dwBufferLength = ARRAY_SIZE * sizeof(waveIn[0]);
    WaveInHdr.dwBytesRecorded = 0;
    WaveInHdr.dwUser = 0L;
    WaveInHdr.dwFlags = 0L;
    WaveInHdr.dwLoops = 0L;

    result = waveInPrepareHeader(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));
    if (result != MMSYSERR_NOERROR) {
        cout << "Error preparing header." << endl;
        waveInClose(hWaveIn);
        return;
    }

    result = waveInAddBuffer(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));
    if (result != MMSYSERR_NOERROR) {
        cout << "Error adding buffer." << endl;
        waveInUnprepareHeader(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));
        waveInClose(hWaveIn);
        return;
    }

    // Start recording
    result = waveInStart(hWaveIn);
    if (result != MMSYSERR_NOERROR) {
        cout << "Error starting recording." << endl;
        waveInUnprepareHeader(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));
        waveInClose(hWaveIn);
        return;
    }

    cout << "Recording for " << RECORD_SECONDS << " seconds..." << endl;

    // Wait until recording completes (poll the buffer flags)
    while (!(WaveInHdr.dwFlags & WHDR_DONE)) {
        Sleep(10);  // Wait and check every 10 ms
    }

    // Stop recording
    waveInStop(hWaveIn);
    waveInUnprepareHeader(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));
    waveInClose(hWaveIn);

    cout << "Recording complete." << endl;
}

void LiveRecording::SaveToFile() {
    std::ofstream outFile("live_amps.txt");
    if (!outFile.is_open()) {
        std::cout << "Error opening file for writing." << std::endl;
        return;
    }

    for (int i = 0; i < ARRAY_SIZE; ++i) {
        outFile << waveIn[i] << "\n";
    }
    outFile.close();
    std::cout << "Amplitudes saved to live_amps.txt" << std::endl;
}

std::vector<double> LiveRecording::readLiveData() {
    std::vector<double> data;
    std::ifstream file("live_amps.txt");
    double value;

    while (file >> value) {
        data.push_back(value);
    }
    return data;
}

void LiveRecording::normalizeData(std::vector<double>& data) {
    double sum = 0.0;
    for (size_t i = 0; i < data.size(); ++i) {
        sum += data[i];
    }
    double mean = sum / data.size();

    double maxAmp = 0.0;
    for (size_t i = 0; i < data.size(); ++i) {
        data[i] -= mean;
        if (fabs(data[i]) > maxAmp) {
            maxAmp = fabs(data[i]);
        }
    }

    if (maxAmp > 0) {
        for (size_t i = 0; i < data.size(); ++i) {
            data[i] /= maxAmp;
        }
    }
}

void LiveRecording::applyHammingWindow(std::vector<double>& frame) {
    for (size_t i = 0; i < frame.size(); ++i) {
        frame[i] *= 0.54 - 0.46 * cos(2 * M_PI * i / (frame.size() - 1));
    }
}

std::vector<double> LiveRecording::computeAutocorrelation(const std::vector<double>& frame, int p) {
    std::vector<double> R(p + 1, 0.0);
    for (int i = 0; i <= p; ++i) {
        for (size_t j = 0; j < frame.size() - i; ++j) {
            R[i] += frame[j] * frame[j + i];
        }
    }
    return R;
}

std::vector<double> LiveRecording::levinsonDurbin(const std::vector<double>& R) {
    std::vector<double> A(P + 1, 0.0);
    std::vector<std::vector<double>> alpha(P + 1, std::vector<double>(P + 1, 0.0));
    std::vector<double> E(P + 1, 0.0);
    std::vector<double> K(P + 1, 0.0);

    E[0] = R[0];

    for (int i = 1; i <= P; ++i) {
        double sum = 0.0;
        for (int j = 1; j < i; ++j) {
            sum += alpha[i - 1][j] * R[i - j];
        }

        K[i] = (R[i] - sum) / E[i - 1];
        alpha[i][i] = K[i];

        for (int j = 1; j < i; ++j) {
            alpha[i][j] = alpha[i - 1][j] - K[i] * alpha[i - 1][i - j];
        }

        E[i] = (1.0 - K[i] * K[i]) * E[i - 1];
    }

    for (int j = 1; j <= P; ++j) {
        A[j] = alpha[P][j];
    }

    return A;
}

std::vector<double> LiveRecording::computeCepstral(const std::vector<double>& A) {
    std::vector<double> C(P + 1, 0.0);

    for (int m = 1; m <= P; ++m) {
        C[m] = A[m];
        for (int k = 1; k < m; ++k) {
            C[m] += (k * C[k] * A[m - k]) / m;
        }
    }

    return C;
}

void LiveRecording::applyRaisedSineWindow(std::vector<double>& C) {
    for (size_t i = 1; i < C.size(); ++i) {
        C[i] *= (1.0 + (P / 2.0) * sin(M_PI * i / P));
    }
}

double LiveRecording::tokhuraDistance(const std::vector<double>& c1, const std::vector<double>& c2) {
    double distance = 0.0;
    for (size_t i = 0; i < weightsArray.size() && i < c1.size() && i < c2.size(); i++) {
        distance += weightsArray[i] * std::pow(c1[i] - c2[i], 2);
    }
    return distance;
}

void LiveRecording::ProcessRecording() {
    std::vector<Frame> frames;
    std::vector<double> data = readLiveData();
    normalizeData(data);

    for (size_t i = 0; i + FRAME_SIZE <= data.size(); i += HOP_LENGTH) {
        Frame frame;
        frame.data.assign(data.begin() + i, data.begin() + i + FRAME_SIZE);

        applyHammingWindow(frame.data);
        frame.R = computeAutocorrelation(frame.data, P);
        frame.A = levinsonDurbin(frame.R);
        frame.C = computeCepstral(frame.A);
        applyRaisedSineWindow(frame.C);

        frames.push_back(frame);
    }

    std::ofstream outFile("live_cepstral.txt");
    // Replace range-based for loop
    for (size_t i = 0; i < frames.size(); i++) {
        outFile << frames[i].C[1];
        for (size_t j = 2; j <= P; ++j) {
            outFile << "," << frames[i].C[j];
        }
        outFile << "\n";
    }
    outFile.close();
}

int LiveRecording::findClosestCodebookIndex(const std::vector<double>& cepstralCoeffs, 
                                          const std::vector<std::vector<double>>& codebook) {
    double minDistance = (std::numeric_limits<double>::max)();
    int minIndex = 0;

    for (size_t i = 0; i < codebook.size(); i++) {
        double distance = tokhuraDistance(cepstralCoeffs, codebook[i]);
        if (distance < minDistance) {
            minDistance = distance;
            minIndex = static_cast<int>(i);
        }
    }
    return minIndex;
}

void LiveRecording::GenerateObservationSequence() {
    auto codebook = loadCodebook("generatedCodeBook.csv");
    if (codebook.empty()) {
        std::cout << "Error: Could not load codebook" << std::endl;
        return;
    }

    std::ifstream cepstralFile("live_cepstral.txt");
    std::ofstream outputFile("live_observation.txt");
    
    if (!cepstralFile.is_open() || !outputFile.is_open()) {
        std::cout << "Error: Could not open required files" << std::endl;
        return;
    }

    std::string line;
    while (std::getline(cepstralFile, line)) {
        std::stringstream ss(line);
        std::vector<double> coeffs;
        double coeff;
        
        while (ss >> coeff) {
            coeffs.push_back(coeff);
            if (ss.peek() == ',') ss.ignore();
        }
        
        int closestIndex = findClosestCodebookIndex(coeffs, codebook);
        outputFile << closestIndex << " ";
    }
    
    outputFile << std::endl;
    outputFile.close();
    cepstralFile.close();
}