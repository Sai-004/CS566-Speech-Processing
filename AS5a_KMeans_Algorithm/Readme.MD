### Assignment 5a: K-Means Algorithm

## K-Means Algorithm for Speech Vector Quantization

This console application implements the K-Means clustering algorithm using the Tokhura distance measure to classify speech vectors based on cepstral coefficients (Ci). The code processes vector data from a CSV file and clusters them into predefined clusters for vector quantization.

### Usage

1. **Prepare Input Data**: Ensure your input vector data is stored in a `universe.csv` file located in the same directory as the application.
   
2. **Run the Application**: 
   - Build and run the project in Visual Studio.
   - The program will initialize the codebook with random vectors from the universe, perform clustering, and iteratively update the centroids.
   - The final codebook vectors will be printed in the console.

3. **Results**: 
   - The application outputs the average distortion for each iteration.
   - Once convergence is reached, the final codebook vectors are displayed.

4. **Exit**: After displaying the final codebook, the program exits after waiting for any key press.

### Key Functions

- **`ComputeTokhuraDistance()`**: 
  Calculates the Tokhura distance between two vectors using predefined Tokhura weights.

- **`main()`**: 
  Implements the main K-Means clustering loop, initializes the codebook, assigns vectors to clusters, and updates the centroids.

### Notes

1. **Codebook Size**: The application uses 8 clusters by default. Modify the `k` variable to change the number of clusters.
2. **Convergence Threshold**: The convergence is based on a threshold (`DELTA`) of 0.0001. This can be adjusted to refine the stopping criteria.
3. **Random Initialization**: The initial codebook vectors are selected randomly from the input universe. Each run might result in different clusters unless random seeding is fixed.