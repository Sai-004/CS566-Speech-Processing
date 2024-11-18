# CS566 Speech Processing Project: Speech-Activated Music Player Manager (Project Number: 10)

## Project Overview
This application enables music playlist control through voice commands, recognizing five basic commands (create, add, play, stop, shuffle) and digits (0-9) for song selection.

## Directory Structure
```
Project_Root/
├── SpeechProcessingProject/
│   ├── Form1.h                # Main application UI and logic
│   ├── HMMCore.h             # HMM functionality declarations
│   ├── HMMCore.cpp           # HMM functionality implementation
│   ├── LiveRecording.h       # Audio recording declarations
│   ├── LiveRecording.cpp     # Audio recording implementation
│   ├── Test.h               # Test framework declarations
│   ├── Test.cpp             # Test framework implementation
│   ├── Models/              # HMM model files
│   │   ├── command_create.txt
│   │   ├── command_add.txt
│   │   ├── command_play.txt
│   │   ├── command_stop.txt
│   │   ├── command_shuffle.txt
│   │   └── model_[0-9].txt
│   ├── Music/               # Directory for .wav files
│   └── files/               # Initial HMM parameters
│       ├── A_ijMatrix.txt
│       ├── B_jkMatrix.txt
│       └── P_iVector.txt
```

## Setup Instructions
1. Open Visual Studio 2010
2. Load the Main Project 
   - Two projects included in the zip file:
     - `MusicPlaylistManagerTest`: Testing recognition model for commands
     - `SpeechProcessingProject`: Main project with UI and functionality

3. Configure Project Settings:
   - Right-click project → Properties
   - C/C++ → General → Additional Include Directories: Add `$(ProjectDir)` if missing
   - Linker → Input → Additional Dependencies: Add `winmm.lib` if missing

4. Build Solution and Run

## Running the Application
### UI Components
- Left panel: Available songs
- Right panel: Current playlist and playback controls
- Bottom: Test and speech recognition controls

### Voice Commands
- `create` - Creates new playlist
- `add` [number] - Adds song by index
- `play` - Plays selected song
- `stop` - Stops playback
- `shuffle` - Randomizes playlist

### Testing
- Click "Run Tests" button for test suite execution
- Tests verify:
  - Command recognition
  - Digit recognition
  - Complete flow testing

## Important Functions
### Form1.h
- `LoadModels()` - Loads HMM models for commands and digits
- `RecognizeCommand()` - Identifies spoken command
- `RecognizeDigit()` - Identifies spoken digit
- `CreateNewPlaylist()` - Creates new playlist
- `AddSongToPlaylist()` - Adds song to current playlist
- `PlayCurrentSong()` - Plays selected song
- `ShufflePlaylist()` - Randomizes playlist
- `RunTests()` - Executes test codes

### LiveRecording.cpp
- `StartRecording()` - Initiates microphone recording
- `StopRecording()` - Stops recording
- `ProcessRecording()` - Processes audio to observation sequence
- `GenerateObservationSequence()` - Creates HMM observation sequence
- `ListInputDevices()` - Lists audio input devices
- `SaveToFile()` - Saves audio data

### HMMCore.cpp
- `loadModelFromFile()` - Loads HMM parameters
- `evaluateSequence()` - Evaluates observation sequence
- `calculateSequenceProbability()` - Calculates sequence probability
- `loadObservationSequence()` - Loads observations
- `initializeModelFromFiles()` - Initializes HMM model
- `averageModels()` - Averages multiple HMM models

## Known Limitations
1. WAV files only supported
2. Background noise affects recognition accuracy
3. Single active playlist

## Testing Instructions
1. Use "Run Tests" button
2. Check TestData folder for generated data
3. Monitor console for test progress
4. Review test results in message box