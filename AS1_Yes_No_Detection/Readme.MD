## Assignment 1: Yes or No Detection

This console application performs speech detection to classify audio as `yes` or `no` based on *Zero Crossing Rate* (ZCR) and *energy thresholds*. 
It provides options to determine the ZCR threshold from samples or to test a specific audio file to decide `yes` or `no`.

## Usage:
On running in the visual studio, console application shows with the options to the user.
1. **Process sample audio files and get ZCR Threshold:** Select this option to use the sample files in the `audio` folder for calculating the ZCR and energy thresholds based on the data and save the results in `results/training_results.csv`.
2. **Test an Audio and Decide 'Yes' or 'No':** Select this option to test an audio file from the `test` folder. The application will use the calculated ZCR threshold to classify the audio as `yes` or `no` and save the result in `results/test_output.csv`.
3. **Exit:** Exits the application.

## Functions Used:
- `calculateEnergy(short* frame, int frameSize)`  
  Calculates the energy of a given frame of audio files. Energy is computed as the sum of squared sample values.

- `calculateZCR(short* frame, int frameSize)`
  Calculates the Zero Crossing Rate (ZCR) of a given frame of audio files. ZCR is the rate at which the signal changes sign and in the function its calculated per sample in a frame.

- `loadAudioData(const char* fileName, short* audioData)`
  Loads audio data from a text file into a static array. Skips the first 5 lines of metadata if present and reads audio sample values.

- `train()`
  Processes sample audio files to calculate average ZCR and energy values for `yes` and `no` categories. Determines the ZCR and energy thresholds and saves the results to `results/training_results.csv`.

- `test(double zcrThreshold)`
  Tests an audio file against the trained ZCR threshold. Computes the average ZCR of the test file and makes a decision based on the threshold. Saves the test results to `results/test_output.csv`.
  Make sure to replace the file name with required audio file name to test it.

## Folder Structure:

**audio:** Contains sample audio files in text format (e.g., `yes_1.txt`, `no_1.txt`). \
**test:** Contains test audio files (e.g., `yes_test.txt`).
results: Stores the output CSV files with the results (e.g., `training_results.csv`, `test_output.csv`).

## Main Menu Options:
1. **Process sample audio files and get ZCR Threshold:** Process samples and updates the ZCR threshold.
2. **Test an Audio and Decide 'Yes' or 'No':** Tests a specified audio file and provides a decision based on the calculated ZCR threshold.
3. **Exit:** Exits the application.

## Example Usage:
1. Run the application.
2. Select option 1 to process samples and get the ZCR threshold.
3. After threshold is found, select option 2 to test an audio file.
4. Decision will be displayed on console.
5. View results in the `results` folder.

---

## Note: 
1. If need to test on new files, trim audio to the word and place the .txt format in the test folder. Make sure to replace the name of file in the code in `test` function to test other files.
2. You will find that option 2 wont work unless you use option 1 atleast once. I have implemented it to show the process of finding thresholds first and then test which is the actual process. 
   *This can be changed if not required by just removing the if case in the option 2 and replacing by calculated threshold.*
3. ZCR values calculated are per sample (divided by frameSize) so they are smaller values (e.g. 0.19, 0.06).
4. When converted to `.txt` in *Cooledit 2000*, metadata is present in the file. It is handled in the code. If no metadata in the new testing files, comment out the skipping lines part in the code (commented in code).