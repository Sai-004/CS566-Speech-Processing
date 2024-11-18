### Assignment 6: Hidden Markov Models (HMM) – Baum-Welch Algorithm for Parameter Re-estimation

This console assignment implements the Baum-Welch algorithm for re-estimating HMM parameters and calculating the probability of the observation sequence at each iteration.

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Problem 3: Baum-Welch Algorithm for Parameter Re-estimation

## Usage:

1. **Input Data**: 
   - The matrix files `A_ijMatrix.txt`, `B_jkMatrix.txt`, and `P_iVector.txt` should be located in the `files` folder in the same directory as the application.
   - The observation sequence should be provided in `ObservationSequence.txt` in the same folder `test`.

2. **Run the Application**:
   - Open the `210101070_HMM3.sln` solution in Visual Studio 2010 from `210101070_HMM3` folder and locate the `.cpp` file.
   - Build and run the project in Visual Studio.

3. **Results**:
   - The program outputs the initial matrices and vectors.
   - For each iteration, it displays:
     * The sequence probability
     * The probability improvement from the previous iteration
     * Updated A matrix, B matrix, and Pi vector
   - The algorithm continues until convergence or maximum iterations are reached.

## Key Functions:

- `baumWelch()`: Implements the Baum-Welch algorithm for re-estimating HMM parameters.
- `calculateSequenceProbability()`: Computes the probability of the observation sequence at each iteration.
- `loadMatrixFromFile()` and `loadArrayFromFile()`: Load matrix and array data from input files.
- `displayMatrix()`: Displays matrices for debugging and output purposes.

## Observations:

- The program shows probability improvement at each iteration.
- Convergence is reached when probability improvement becomes zero.
- Re-estimated parameters (A, B, Pi) are displayed after each iteration.
- The algorithm implements both forward and backward procedures as part of the Baum-Welch algorithm.

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## Notes:
1. Matrix and Vector Sizes:
   - `A` matrix size: 5 x 5
   - `B` matrix size: 5 x 32
   - `Pi` vector size: 5

2. Implementation Details:
   - Uses dynamic arrays instead of vectors for better debugging compatibility with Visual Studio 2010
   - Includes proper memory management with allocation and deallocation
   - The algorithm stops when probability improvement becomes zero. Maximum iterations can be specified (default: 100)

3. Output Format:
   - Matrices are displayed with 4 decimal place precision
   - Probabilities are shown in scientific notation with 10 decimal places
   - Probability improvements are calculated and displayed for each iteration