### Assignment 5b: LBG Algorithm

## Linde-Buzo-Gray (LBG) Algorithm for Vector Quantization

This console application implements the Linde-Buzo-Gray (LBG) algorithm for vector quantization using cepstral coefficients (Ci). The algorithm involves an iterative process of splitting clusters and applying K-Means clustering to generate an optimal codebook for vector quantization.

### Usage

1. **Prepare Input Data**: Ensure your input vector data is stored in a `universe.csv` file located in the same directory as the application.

2. **Run the Application**: 
   - Build and run the project in Visual Studio.
   - The program first calculates the initial centroid of the universe, then iteratively splits the codebook vectors and applies K-Means clustering until the desired codebook size is reached.
   - The codebook vectors are displayed after each split and K-Means iteration.

3. **Results**: 
   - The application outputs the codebook size after each split and the average distortion for each K-Means iteration.
   - The final codebook vectors are displayed at the end.

4. **Exit**: After processing, the program exits after waiting for any key press.

### Key Functions

- **`ComputeTokhuraDistance()`**: 
  Computes the Tokhura distance between two vectors based on predefined Tokhura weights.

- **`kMeans()`**: 
  Runs the K-Means clustering algorithm on a set of vectors and updates the codebook iteratively.

- **`main()`**: 
  Implements the LBG algorithm by splitting clusters and running K-Means clustering until the desired codebook size is achieved.

### Notes

1. **Codebook Size**: The desired codebook size is set to 8 by default. This can be modified by changing the `k` variable in the code.
2. **Split Factor**: The splitting factor (`EPSILON`) is set to 0.03. Adjust this value if different splitting behavior is desired.
3. **Initial Centroid**: The LBG algorithm starts by calculating the centroid of the entire universe as the initial codebook vector.
4. **Convergence Threshold**: K-Means convergence is based on a threshold (`DELTA`) of 0.0001, which can be adjusted for finer control of the iterations.