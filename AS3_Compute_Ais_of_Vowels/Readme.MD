## Assignment 3: Calculate ai's from a speech signal

This console application performs *Linear Predictive Coding* (LPC) analysis on speech signals. 
It calculates the LPC coefficients (*`Ai`*) and autocorrelation coefficients (*`Ri`*) for the 5 continuous steady frames with maximum energy in the input signal.

## Usage:

1. **Prepare Input:** Ensure you have a text file containing the speech signal samples (one sample per line).

2. **Run the Application:** Execute the program. When prompted, enter the path to your input file.

3. **View Results:** The program processes the audio and writes the results to two files in the "Results" folder:
   - `Ai_values.txt`: Contains the LPC coefficients for each selected frame.
   - `Ri_values.txt`: Contains the autocorrelation coefficients for each selected frame.

4. **Exit:** After processing, the program displays a completion message and exits automatically.

## Key Functions:

- `ComputeRi(double* frame, int frameSize, double* Ri)`  
  Computes the autocorrelation coefficients (*`Ri`*) for a given frame of audio data.

- `ComputeAi(double* R, double* A)`  
  Calculates the LPC coefficients (*`Ai`*) using the Levinson-Durbin algorithm based on the autocorrelation coefficients.

- `applyNormalization(double* data, int dataSize)`  
  Performs DC shift and amplitude normalization on the entire input signal.

- `applyHammingWindow(double* frame, int frameSize)`  
  Applies a Hamming window to a frame of audio data to reduce spectral leakage.

- `calculateFrameEnergy(double* frame, int frameSize)`  
  Calculates the energy of a frame of audio data.

## Process Overview:

1. **Data Loading:** The program reads the input signal from the specified file.

2. **Normalization:** The entire signal undergoes DC shift and amplitude normalization.

3. **Energy Calculation:** The energy for each frame in the signal is computed.

4. **Frame Selection:** The program identifies the frame with maximum energy and selects 5 continuous frames around it (2 before, the max energy frame, and 2 after).

5. **LPC Analysis:** For each selected frame:
   - A Hamming window is applied.
   - Autocorrelation coefficients (*`Ri`*) are computed.
   - LPC coefficients (*`Ai`*) are calculated using the Levinson-Durbin algorithm.

6. **Result Output:** The *`Ri`* and *`Ai`* values for each selected frame are written to their respective output files.
---
## Notes:
1. **Input Format:** The input file should contain one sample per line, with samples represented as floating-point numbers.
2. **Frame Size:** The program uses a frame size of 320 samples. Adjust the `FRAME_SIZE` constant if a different frame size is required.
3. **LPC Order:** The current implementation uses an LPC order of 12. Modify the `p` constant to change this.
4. **Output:** The program outputs both *`Ai`* and *`Ri`* values for each of the 5 selected frames, along with their frame numbers and offsets from the max energy frame.
5. **Error Handling:** The program includes basic error checking for file operations. Ensure you have write permissions in the "Results" folder.