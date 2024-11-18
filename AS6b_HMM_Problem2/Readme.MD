### Assignment 6: Hidden Markov Models (HMM) – Forward-Backward Algorithm and Viterbi Algorithm

This console assignment comprises of implementing two key algorithms in HMM: the Forward-Backward algorithm for evaluating the probability of an observation sequence and the Viterbi algorithm for finding the most likely state sequence.

The code processes matrix data from text files to perform these operations based on predefined states and observations.

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Problem 2: Viterbi Algorithm for State Sequence Prediction

## Usage:

1. **Input Data**:
   - The same input files (`A_ijMatrix.txt`, `B_jkMatrix.txt`, `P_iVector.txt`, and `ObservationSequence.txt`) are used.

2. **Run the Application**:
    - Open the `210101070_HMM2.sln` solution in Visual Studio 2010 from `210101070_HMM2` folder and locate the `.cpp` file.
    - Build and run the project in Visual Studio.

3. **Results**:
   - The program outputs the matrices, vectors, and the observation sequence.
   - It calculates the most likely state sequence using the Viterbi algorithm and prints the sequence in both 0-based and 1-based indexing formats.

## Key Functions:

- **`ViterbiAlgorithm()`**: Implements the Viterbi algorithm for finding the most probable sequence of hidden states based on the observation sequence.
- **`loadMatrixFromFile()`** and **`loadVectorFromFile()`**: Load matrix and vector data from input files.

## Observations:

- Length of the state sequence: 85
- Most likely state sequence: (0 based indexing):
0 1 2 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4

- Most likely state sequence: (1 based indexing):
1 2 3 4 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Notes:
1. Matrix and Vector Sizes:
   - `A` matrix size: 5 x 5
   - `B` matrix size: 5 x 32
   - `Pi` vector size: 5
   - Observation sequence length: 85

2. Consistency: The Forward and Backward algorithms in Problem 1 produce identical results for the observation probability (`P(O|λ)`).

3. State Sequence: The Viterbi algorithm in Problem 2 outputs both 0-based and 1-based indexed state sequences for better interpretation.

