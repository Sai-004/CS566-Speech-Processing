#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Windows.h>

#define FRAME_SIZE 400
#define ARRAY_SIZE 16025 * 3
#define ENERGY_THRESHOLD 10000

short int waveIn[ARRAY_SIZE];

#pragma comment(lib, "winmm.lib")

void StartRecord()
{
    const int NUMPTS = 16025 * 3;
    int sampleRate = 16025;
    HWAVEIN hWaveIn;
    MMRESULT result;
    WAVEFORMATEX pFormat;
    pFormat.wFormatTag = WAVE_FORMAT_PCM;
    pFormat.nChannels = 1;
    pFormat.nSamplesPerSec = sampleRate;
    pFormat.nAvgBytesPerSec = sampleRate * 2;
    pFormat.nBlockAlign = 2;
    pFormat.wBitsPerSample = 16;
    pFormat.cbSize = 0;

    result = waveInOpen(&hWaveIn, WAVE_MAPPER, &pFormat, 0L, 0L, WAVE_FORMAT_DIRECT);
    WAVEHDR WaveInHdr;
    WaveInHdr.lpData = (LPSTR)waveIn;
    WaveInHdr.dwBufferLength = NUMPTS * 2;
    WaveInHdr.dwBytesRecorded = 0;
    WaveInHdr.dwUser = 0L;
    WaveInHdr.dwFlags = 0L;
    WaveInHdr.dwLoops = 0L;

    waveInPrepareHeader(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));
    result = waveInAddBuffer(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));
    result = waveInStart(hWaveIn);

    printf("Recording started... (3 seconds)\n");
    Sleep(3000);

    waveInClose(hWaveIn);
    printf("Recording Saved.\n");
}

double calculateZCR(const int *frame, int start, int end)
{
    int zcr_count = 0;
    for (int i = start + 1; i <= end; ++i)
    {
        if ((frame[i - 1] > 0 && frame[i] < 0) || (frame[i - 1] < 0 && frame[i] > 0))
        {
            zcr_count++;
        }
    }
    return (double)zcr_count / (end - start + 1);
}

double calculateEnergy(const int *frame, int start)
{
    double energy = 0;
    for (int i = start; i < start + FRAME_SIZE; ++i)
    {
        energy += (double)frame[i] * frame[i];
    }
    printf("Frame starting at %d has energy: %.2f\n", start, energy);
    return energy;
}

void processAudio()
{
    int frame[ARRAY_SIZE];
    for (int i = 0; i < ARRAY_SIZE; ++i)
        frame[i] = static_cast<int>(waveIn[i]);

    int size = ARRAY_SIZE;
    int start = 0, end = size - 1, startTrue = 0;
    for (int i = 0; i < size - FRAME_SIZE; i += FRAME_SIZE)
    {
        double energy = calculateEnergy(frame, i);

        if (!startTrue && energy > ENERGY_THRESHOLD)
        {
            start = i;
            startTrue = 1;
            printf("Speech starts at frame: %d\n", start);
        }

        if (startTrue && energy < ENERGY_THRESHOLD)
        {
            end = i;
            printf("Speech ends at frame: %d\n", end);
            break;
        }
    }

    double zcr = calculateZCR(frame, start, end);
    printf("ZCR detected: %.5f\n", zcr);
    printf("Word detected: %s\n", (zcr > 0.07) ? "Yes" : "No");
}

int main()
{
    StartRecord();
    processAudio();
    printf("Exiting program.\n");
    Sleep(1000);
    return 0;
}