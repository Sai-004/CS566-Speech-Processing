### Assignment 8: Hidden Markov Model (HMM) Digit Recognition System
A comprehensive implementation of speech-based digit recognition using Hidden Markov Models with Baum-Welch algorithm for parameter re-estimation and live recognition capabilities.

## Overview

This system implements speech recognition for digits 0-9 using:
- Hidden Markov Models (HMM) for pattern recognition
- Baum-Welch algorithm for parameter re-estimation
- Live audio recording and processing functionality
- Real-time digit recognition capabilities


## System Architecture

1. **HMM Implementation**
- State Count: 5 states per model
- Maximum Sequence Length: 100 observations
- Observation Symbols: 32 distinct symbols
- Training Set: 30 samples per digit
- Test Set: 10 samples per digit

2. **Live Recording Module**
- Sample Rate: 16025 Hz
- Frame Size: 320 samples
- Hop Length: 240 samples
- LPC Order: 12 coefficients
- Recording Duration: 1.2 seconds

## Key Features
### Model Training and Testing
1. **Training Phase**
   - Uses Baum-Welch algorithm for parameter re-estimation
   - 60 iterations per training sequence
   - Model averaging across training samples
   - Automatic model saving for future use

2. **Testing Phase**
   - Cross-validation against test dataset
   - Uncertainty detection based on probability differences
   - Detailed accuracy metrics per digit
   - Comprehensive error analysis

### Live Recognition
1. **Audio Capture**
   - Real-time audio recording
   - Device selection capability
   - Automatic gain control

2. **Signal Processing**
   - Hamming window application
   - Autocorrelation computation
   - Levinson-Durbin recursion
   - Cepstral coefficient extraction

## Usage
1. **Initial Setup**
Required folder structure
files/
  ├── A_ijMatrix.txt    Initial transition matrix (5x5)
  ├── B_jkMatrix.txt    Initial emission matrix (5x32)
  └── P_iVector.txt     Initial state probabilities (5)

models/                 Trained model storage
210101070_dataset/     Training and test data

2. **Running the System**
The system provides a menu-driven interface with the following options:
1. Train New Models
2. Load Existing Models
3. Test Models
4. Live Recognition
5. Exit

3. **Model Training**
- Loads initial parameters from files/
- Iteratively improves using Baum-Welch algorithm
- Automatically saves trained models to models/

4. **Live Recognition**
- Speak a single digit
- System processes and displays recognition result

## Implementation Details
### Core Functions
1. **Training Functions**
   - `trainModels()`: Main training orchestration
   - `baumWelch()`: Parameter re-estimation
   - `averageModels()`: Model combination

2. **Testing Functions**
   - `testModels()`: Evaluation on test set
   - `evaluateSequence()`: Single sequence probability calculation

3. **Live Processing**
   - `StartRecording()`: Audio capture
   - `ProcessRecording()`: Feature extraction
   - `GenerateObservationSequence()`: Symbol sequence generation

## Technical Specifications
- Matrix Dimensions:
  - Transition Matrix (A): 5×5
  - Emission Matrix (B): 5×32
  - Initial Probabilities (π): 5×1
- Minimum Probability: 1e-15
- Convergence threshold: Based on probability improvement

---------------------------------------------------------------------------------------------------------------------------------------------------------------
### Notes
1. The system uses dynamic memory allocation for matrix operations
2. Implements proper error handling for file operations and memory management
3. Includes uncertainty detection for more reliable recognition