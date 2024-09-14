## Assignment 2: Yes or No Detection 2

This console application performs real-time speech detection to classify audio as `yes` or `no` based on *Zero Crossing Rate* (ZCR) and *energy thresholds*. 
The user can record an audio sample directly within the console, and the application will process the recorded audio to make a classification.

## Usage:
1. **Run the Application:** The application will begin by recording a 3-second audio sample.
2. **View Results:** The program processes the audio and displays whether the word spoken was `yes` or `no` based on the calculated ZCR and energy.
3. **Exit:** After displaying the results, the program exits automatically.

## Functions Used:
- `StartRecord()`  
  This function records 3 seconds of audio at a sample rate of 16 kHz using the system’s microphone. The recorded audio data is stored in a buffer for further processing.

- `calculateEnergy(const int* frame, int start)`  
  Calculates the energy of a frame of audio data. Energy is computed as the sum of squared sample values for the given frame.

- `calculateZCR(const int* frame, int start, int end)`  
  Calculates the Zero Crossing Rate (ZCR) of a segment of the audio data. ZCR is determined by counting the number of times the signal changes sign. (from +ve to -ve or -ve to +ve)

- `processAudio()`  
  Processes the recorded audio by calculating the energy and ZCR. It determines the start and end of the speech in the recorded audio based on the energy threshold. 
  Then, it calculates the ZCR for the identified speech segment to classify the word as `yes` or `no.`

---
## Notes:
1. **Recording:** The application now records audio directly within the console, eliminating the need to use pre-recorded samples.
2. **Thresholds:** The ZCR and energy thresholds are used to detect and classify the speech. The current thresholds assume typical values for `yes` or `no.`
3. **Output:** The application prints the start and end frames of the detected speech, the calculated ZCR, and the classification result (`Yes` or `No`).