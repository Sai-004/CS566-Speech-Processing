# Assignment 4: Calculate Ai's from a Speech Signal

This console application performs **vowel recognition** using **Linear Predictive Coding (LPC)** and **Tokhura distance** based on the speech signal's cepstral coefficients. The program processes speech samples, calculates LPC coefficients (Ai), autocorrelation coefficients (Ri), and cepstral coefficients (Ci), then predicts the vowel based on Tokhura distance between test and training data.

---

## Features
- **LPC Analysis**: Calculates LPC coefficients (Ai) using Levinson-Durbin Algorithm.
- **Cepstral Coefficients (Ci)**: Derived from Ai and smoothed using a sine window.
- **Tokhura Distance**: Used for vowel prediction by comparing Ci's of test data to reference values from the training set.
- **Vowel Recognition**: Supports vowel recognition for 'a', 'e', 'i', 'o', 'u' using preprocessed data from training files.

---

## Usage

### 1. Prepare Input:
Ensure that you have the following files in the appropriate directories:
- **Training data**: Text files in the `TrainData` folder, where each file contains speech signal samples (one sample per line) for different vowels (e.g., `210101070_a_1.txt`).
- **Test data**: Text files in the `TestData` folder for validation and vowel prediction (e.g., `210101070_a_21.txt`).

### 2. Run the Application:
Execute the program. The application will automatically:
- Process the training files.
- Generate reference Ci values for each vowel.
- Predict the vowels from test files using Tokhura distance.

### 3. View Results:
- The program will output the predicted vowel for each test file and display Tokhura distance values for each vowel ('a', 'e', 'i', 'o', 'u').
- It will also print the accuracy for each vowel and the overall accuracy.

### 4. Exit:
The program will display a completion message and automatically exit.

---

## Key Functions

### `ComputeRi()`
- Computes the autocorrelation coefficients (Ri) for each frame in the audio signal using stable frames of maximum energy.

### `ComputeAi()`
- Calculates the Linear Predictive Coefficients (Ai) using the Levinson-Durbin algorithm based on the autocorrelation coefficients (Ri).

### `ComputeCi()`
- Computes the cepstral coefficients (Ci) based on the LPC coefficients (Ai), applying a sine window for smoothing.

### `applyHammingWindow()`
- Applies a Hamming window to the stable frames of audio data to reduce spectral leakage.

### `saveAvgCi()`
- Saves the average cepstral coefficients for each vowel to reference files for later comparison during vowel prediction.

### `calculateTokhuraDistance()`
- Calculates the Tokhura distance between the test file's Ci's and reference Ci's for each vowel.

---

## Process Overview

1. **Data Loading**: 
   - The program reads the speech samples from the training and test files.
   
2. **Normalization**:
   - Each signal undergoes **DC shift** and **amplitude normalization**.

3. **Frame Selection**:
   - The program computes the energy of each frame in the signal and identifies 5 continuous **stable frames** around the frame with maximum energy.

4. **LPC and Cepstral Analysis**:
   - For each of the 5 stable frames:
     - A **Hamming window** is applied to reduce spectral leakage.
     - **Autocorrelation coefficients (Ri)** are computed.
     - **LPC coefficients (Ai)** are calculated using the Levinson-Durbin algorithm.
     - **Cepstral coefficients (Ci)** are derived from the LPC coefficients and smoothed using a **sine window**.

5. **Training and Reference Generation**:
   - The program processes the training files for each vowel ('a', 'e', 'i', 'o', 'u') and saves the **average Ci** values in reference files.

6. **Vowel Prediction**:
   - The program processes each test file and calculates the **Tokhura distance** between the test Ci values and the reference Ci values for each vowel.
   - The vowel with the **minimum Tokhura distance** is predicted for the test file.

7. **Result Output**:
   - The program displays the Tokhura distance for each vowel and prints the predicted vowel.
   - It also calculates the **accuracy** for each vowel and prints the **overall accuracy**.

---

## Notes
1. **Input Format**: The input files should contain one sample per line, represented as integers or floating-point numbers.
2. **Training Data**: The program processes **20 training files per vowel**.
3. **Testing**: Each vowel has **10 test files** for evaluation.
4. **Output**: 
   - The program outputs both Tokhura distances and predicted vowels, along with the accuracy per vowel and overall.
   - Reference Ci files are stored in the `Reference` folder.
5. **Error Handling**: Basic error checking is included for file operations. Ensure the `TrainData`, `TestData`, and `Reference` directories are properly structured and writable.

---

## Example

### Training Data File Structure
```
TrainData/
|-- 210101070_a_1.txt
|-- 210101070_a_2.txt
...
|-- 210101070_u_20.txt
```

### Test Data File Structure
```
TestData/
|-- 210101070_a_21.txt
|-- 210101070_a_22.txt
...
|-- 210101070_u_30.txt
```

### Reference File Output
```
Reference/
|-- reference_ci_a.txt
|-- reference_ci_e.txt
...
|-- reference_ci_u.txt
```

---

## Adjustments
- **Frame Size**: You can modify the constant `FRAME_SIZE` (currently 320) if a different frame size is needed.
- **LPC Order**: The program uses an LPC order of 12. You can change the `P` constant to adjust the number of prediction coefficients.
- **Test and Training Data Size**: The training data size is set to 20 files per vowel, and testing uses 10 files per vowel. You can adjust these values by modifying the `TRAINDATA_SIZE` constant and the loop logic in `train()` and `test()`.

---


