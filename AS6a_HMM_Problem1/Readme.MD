### Assignment 6: Hidden Markov Models (HMM) – Forward-Backward Algorithm and Viterbi Algorithm

This console assignment comprises of implementing two key algorithms in HMM: the Forward-Backward algorithm for evaluating the probability of an observation sequence and the Viterbi algorithm for finding the most likely state sequence.

The code processes matrix data from text files to perform these operations based on predefined states and observations.

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Problem 1: Forward-Backward Algorithm for Probability Computation

## Usage:

1. **Input Data**: 
   - The matrix files `A_ijMatrix.txt`, `B_jkMatrix.txt`, and `P_iVector.txt` should be located in the `files` folder in the same directory as the application.
   - The observation sequence should be provided in `ObservationSequence.txt` in the same folder `test`.

2. **Run the Application**:
   - Open the `210101070_HMM1.sln` solution in Visual Studio 2010 from `210101070_HMM1` folder and locate the `.cpp` file.
   - Build and run the project in Visual Studio.

3. **Results**:
   - The program outputs the matrices and vectors used.
   - It calculates the probability of the observation sequence using both the Forward and Backward algorithms, with both methods producing the same result.

## Key Functions:

- **`computeForwardProbability()`**: Implements the forward procedure for computing the probability of the observation sequence.
- **`computeBackwardProbability()`**: Implements the backward procedure for computing the same probability as a cross-check.
- **`loadMatrixFromFile()`** and **`loadVectorFromFile()`**: Load matrix and vector data from input files.

## Observations:

- P(O|λ) calculated in both the forward and backward process.
- Forward algorithm probability: 1.15412e-128
- Backward algorithm probability: 1.15412e-128
- Both methods produce the same result.

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Notes:
1. Matrix and Vector Sizes:
   - `A` matrix size: 5 x 5
   - `B` matrix size: 5 x 32
   - `Pi` vector size: 5
   - Observation sequence length: 85

2. Consistency: The Forward and Backward algorithms in Problem 1 produce identical results for the observation probability (`P(O|λ)`).

3. State Sequence: The Viterbi algorithm in Problem 2 outputs both 0-based and 1-based indexed state sequences for better interpretation.

