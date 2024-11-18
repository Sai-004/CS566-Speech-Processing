#pragma once

#include "stdafx.h"
#include "HMMCore.h"
#include "Test.h"
#include "LiveRecording.h"
#include <string>
#include <vcclr.h>

namespace SpeechProcessingProject {

    using namespace System;
    using namespace System::ComponentModel;
    using namespace System::Collections;
    using namespace System::Windows::Forms;
    using namespace System::Data;
    using namespace System::Drawing;
    using namespace System::IO;
    using namespace System::Media;

    public ref class Form1 : public System::Windows::Forms::Form
    {
    private:
        // Native wrapper class to handle unmanaged resources
        ref class NativeResources
        {
        public:
            NativeResources() : commandModels(nullptr), digitModels(nullptr)
            {
                commandModels = new HMMParameters[5];
                digitModels = new HMMParameters[10];
            }

            ~NativeResources()
            {
                this->!NativeResources();
            }

            !NativeResources()
            {
                if (commandModels != nullptr)
                {
                    delete[] commandModels;
                    commandModels = nullptr;
                }
                if (digitModels != nullptr)
                {
                    delete[] digitModels;
                    digitModels = nullptr;
                }
            }

            HMMParameters* commandModels;
            HMMParameters* digitModels;
        };

    public:
        Form1(void)
        {
            InitializeComponent();
            APPLICATION_FOLDER = Application::StartupPath;
            MUSIC_FOLDER = Path::Combine(APPLICATION_FOLDER, "Music");
            MODELS_FOLDER = Path::Combine(APPLICATION_FOLDER, "Models");

            isListening = false;
            recorder = new LiveRecording();
            playlists = gcnew System::Collections::Generic::List<String^>();
            nativeRes = gcnew NativeResources();
            
            // Load HMM models
            LoadModels();
            // Load songs from music folder
            LoadSongs();
        }

		// Method to simulate command processing for testing
		void SimulateCommand()
		{
			RecognitionTimer_Tick(nullptr, nullptr);
		}

		// Method to run all tests
		void RunTests()
		{
			String^ testDir = Path::Combine(APPLICATION_FOLDER, "TestData");
			if (!Directory::Exists(testDir))
			{
				Directory::CreateDirectory(testDir);
			}

			std::string testDataDir = ConvertToStdString(testDir);

			// Generate test data
			if (!Test::GenerateTestObservations(testDataDir))
			{
				MessageBox::Show("Failed to generate test data", "Test Error");
				return;
			}

			// Run individual tests
			bool commandTestResult = Test::RunCommandRecognitionTest(
				nativeRes->commandModels, testDataDir);
        
			bool digitTestResult = Test::RunDigitRecognitionTest(
				nativeRes->digitModels, testDataDir);
        
			bool flowTestResult = Test::RunCompleteFlowTest(
				this, testDataDir);

			// Show test results
			String^ results = String::Format(
				"Test Results:\n\n" +
				"Command Recognition: {0}\n" +
				"Digit Recognition: {1}\n" +
				"Complete Flow: {2}",
				commandTestResult ? "PASSED" : "FAILED",
				digitTestResult ? "PASSED" : "FAILED",
				flowTestResult ? "PASSED" : "FAILED"
			);

			MessageBox::Show(results, "Test Results");
		}

    protected:
        ~Form1()
        {
            if (components)
            {
                delete components;
            }
            if (recorder != nullptr)
            {
                delete recorder;
            }
            if (nativeRes != nullptr)
            {
                delete nativeRes;
            }
        }

    private:
        // UI Controls
        System::ComponentModel::Container ^components;
        System::Windows::Forms::ListBox^ songListBox;
        System::Windows::Forms::ListBox^ playlistListBox;
        System::Windows::Forms::Label^ currentPlaylistLabel;
        System::Windows::Forms::Label^ nowPlayingLabel;
        System::Windows::Forms::Button^ startListeningButton;
        System::Windows::Forms::Timer^ recognitionTimer;
        System::Windows::Forms::Label^ statusLabel;
        System::Windows::Forms::Label^ commandLabel;
        System::Media::SoundPlayer^ mediaPlayer;
		System::Windows::Forms::Button^ testButton;

        // Constants for file paths
        String^ APPLICATION_FOLDER;
        String^ MUSIC_FOLDER;
        String^ MODELS_FOLDER;
        
        // State variables
        bool isListening;
        LiveRecording* recorder;
        NativeResources^ nativeRes;
        String^ currentPlaylist;
        System::Collections::Generic::List<String^>^ playlists;

        // Helper methods for string conversion
		private:
			String^ ConvertToSystemString(const std::string& str) {
				return gcnew String(str.c_str());
			}

			String^ ConvertToSystemString(const char* str) {
				return gcnew String(str);
			}

			std::string ConvertToStdString(String^ str) {
				if (str == nullptr)
					return std::string();
            
				pin_ptr<const wchar_t> wch = PtrToStringChars(str);
				size_t origsize = wcslen(wch) + 1;
				size_t convertedChars = 0;
				const size_t newsize = origsize * 2;
				char* nstring = new char[newsize];
				wcstombs_s(&convertedChars, nstring, newsize, wch, _TRUNCATE);
				std::string result(nstring);
				delete[] nstring;
				return result;
			}

			void LoadModels()
			{
				try
				{
					std::string modelBasePath = ConvertToStdString(MODELS_FOLDER);

					// Verify directory exists
					if (!Directory::Exists(MODELS_FOLDER))
					{
						MessageBox::Show("Models directory not found: " + MODELS_FOLDER, 
							"Error", MessageBoxButtons::OK, MessageBoxIcon::Error);
						return;
					}

					// Load command models
					const char* commands[] = {"create", "add", "play", "stop", "shuffle"};
					for (int i = 0; i < 5; i++)
					{
						std::string command(commands[i]); // Convert const char* to std::string
						std::string fullPath = modelBasePath + "\\command_" + command + ".txt";

						// Verify file exists
						if (!File::Exists(gcnew String(fullPath.c_str())))
						{
							MessageBox::Show("Model file not found: " + gcnew String(fullPath.c_str()),
								"Error", MessageBoxButtons::OK, MessageBoxIcon::Error);
							return;
						}

						if (!loadModelFromFile(nativeRes->commandModels[i], fullPath))
						{
							String^ errorMsg = gcnew String("Failed to load model for command: ");
							errorMsg = String::Concat(errorMsg, gcnew String(command.c_str()));
							MessageBox::Show(errorMsg, "Error", MessageBoxButtons::OK, MessageBoxIcon::Error);
							return;
						}
					}

					// Load digit models
					for (int i = 0; i < 10; i++)
					{
						char digitStr[8];
						_itoa_s(i, digitStr, sizeof(digitStr), 10);
						std::string digit(digitStr); // Convert char array to std::string
						std::string fullPath = modelBasePath + "\\model_" + digit + ".txt";
            
						// Verify file exists
						if (!File::Exists(gcnew String(fullPath.c_str())))
						{
							MessageBox::Show("Model file not found: " + gcnew String(fullPath.c_str()),
								"Error", MessageBoxButtons::OK, MessageBoxIcon::Error);
							return;
						}

						if (!loadModelFromFile(nativeRes->digitModels[i], fullPath))
						{
							MessageBox::Show("Failed to load model for digit: " + i,
								"Error", MessageBoxButtons::OK, MessageBoxIcon::Error);
							return;
						}
					}

					statusLabel->Text = "Status: Models loaded successfully";
					Console::WriteLine("All models loaded successfully");
				}
				catch (Exception^ ex)
				{
					MessageBox::Show("Error loading models: " + ex->Message, "Error",
						MessageBoxButtons::OK, MessageBoxIcon::Error);
				}
			}

			String^ RecognizeCommand()
			{
				const std::string& obsPath = recorder->GetObservationFilePath();
				Console::WriteLine("Evaluating observation sequence from: " + 
					gcnew String(obsPath.c_str()));

				double maxProb = -std::numeric_limits<double>::infinity();
				int bestIndex = 0;
    
				const char* commands[] = {"create", "add", "shuffle", "play", "stop"};
    
				for (int i = 0; i < 5; i++)
				{
					double prob = evaluateSequence(nativeRes->commandModels[i], obsPath);
					Console::WriteLine(String::Format("Probability for {0}: {1}", 
						gcnew String(commands[i]), prob));
        
					if (prob > maxProb)
					{
						maxProb = prob;
						bestIndex = i;
					}
				}
    
				String^ result = gcnew String(commands[bestIndex]);
				Console::WriteLine("Selected command: " + result + " with probability: " + maxProb);
				return result;
			}

			int RecognizeDigit()
			{
				const std::string& obsPath = recorder->GetObservationFilePath();  // Now receives std::string directly
				double maxProb = -std::numeric_limits<double>::infinity();
				int bestDigit = 0;
    
				for (int i = 0; i < 10; i++)
				{
					double prob = evaluateSequence(nativeRes->digitModels[i], obsPath);
					if (prob > maxProb)
					{
						maxProb = prob;
						bestDigit = i;
					}
				}
    
				return bestDigit;
			}

#pragma region Windows Form Designer generated code
        void InitializeComponent(void)
        {
            this->components = gcnew System::ComponentModel::Container();
            
            // Initialize controls
            this->songListBox = (gcnew System::Windows::Forms::ListBox());
            this->playlistListBox = (gcnew System::Windows::Forms::ListBox());
            this->currentPlaylistLabel = (gcnew System::Windows::Forms::Label());
            this->nowPlayingLabel = (gcnew System::Windows::Forms::Label());
            this->startListeningButton = (gcnew System::Windows::Forms::Button());
            this->recognitionTimer = (gcnew System::Windows::Forms::Timer(this->components));
            this->statusLabel = (gcnew System::Windows::Forms::Label());
            this->commandLabel = (gcnew System::Windows::Forms::Label());
            this->mediaPlayer = (gcnew System::Media::SoundPlayer());

            this->SuspendLayout();

            // Form settings
            this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
            this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
            this->ClientSize = System::Drawing::Size(800, 600);
            this->Text = L"Speech Music Player";
            
            // songListBox
            this->songListBox->FormattingEnabled = true;
            this->songListBox->Location = System::Drawing::Point(10, 10);
            this->songListBox->Name = L"songListBox";
            this->songListBox->Size = System::Drawing::Size(380, 540);
            this->Controls->Add(this->songListBox);
            
            // playlistListBox
            this->playlistListBox->FormattingEnabled = true;
            this->playlistListBox->Location = System::Drawing::Point(400, 30);
            this->playlistListBox->Name = L"playlistListBox";
            this->playlistListBox->Size = System::Drawing::Size(380, 360);
            this->Controls->Add(this->playlistListBox);
            
            // currentPlaylistLabel
            this->currentPlaylistLabel->AutoSize = true;
            this->currentPlaylistLabel->Location = System::Drawing::Point(400, 10);
            this->currentPlaylistLabel->Name = L"currentPlaylistLabel";
            this->currentPlaylistLabel->Size = System::Drawing::Size(380, 20);
            this->currentPlaylistLabel->Text = L"Current Playlist: None";
            this->Controls->Add(this->currentPlaylistLabel);
            
            // nowPlayingLabel
            this->nowPlayingLabel->AutoSize = true;
            this->nowPlayingLabel->Location = System::Drawing::Point(400, 400);
            this->nowPlayingLabel->Name = L"nowPlayingLabel";
            this->nowPlayingLabel->Size = System::Drawing::Size(380, 40);
            this->nowPlayingLabel->Text = L"Now Playing: None";
            this->Controls->Add(this->nowPlayingLabel);
            
            // commandLabel
            this->commandLabel->AutoSize = true;
            this->commandLabel->Location = System::Drawing::Point(400, 440);
            this->commandLabel->Name = L"commandLabel";
            this->commandLabel->Size = System::Drawing::Size(380, 30);
            this->commandLabel->Text = L"Last Command: None";
            this->Controls->Add(this->commandLabel);

			// testButton
			this->testButton = (gcnew System::Windows::Forms::Button());
			this->testButton->Location = System::Drawing::Point(10, 560);
			this->testButton->Name = L"testButton";
			this->testButton->Size = System::Drawing::Size(100, 30);
			this->testButton->Text = L"Run Tests";
			this->testButton->Click += gcnew System::EventHandler(this, &Form1::TestButton_Click);
			this->Controls->Add(this->testButton);
            
            // startListeningButton
            this->startListeningButton->Location = System::Drawing::Point(400, 480);
            this->startListeningButton->Name = L"startListeningButton";
            this->startListeningButton->Size = System::Drawing::Size(380, 40);
            this->startListeningButton->Text = L"Start Listening";
            this->startListeningButton->Click += gcnew System::EventHandler(this, &Form1::StartListening_Click);
            this->Controls->Add(this->startListeningButton);
            
            // statusLabel
            this->statusLabel->AutoSize = true;
            this->statusLabel->Location = System::Drawing::Point(400, 530);
            this->statusLabel->Name = L"statusLabel";
            this->statusLabel->Size = System::Drawing::Size(380, 20);
            this->statusLabel->Text = L"Status: Ready";
            this->Controls->Add(this->statusLabel);
            
            // recognitionTimer
            this->recognitionTimer->Interval = 2000;
            this->recognitionTimer->Tick += gcnew System::EventHandler(this, &Form1::RecognitionTimer_Tick);

            this->ResumeLayout(false);
            this->PerformLayout();
        }
#pragma endregion

    private: 
        void LoadSongs()
        {
            try
            {
                EnsureDirectoriesExist();
                array<String^>^ files = Directory::GetFiles(MUSIC_FOLDER, "*.wav");
                songListBox->Items->Clear();
                
                for each (String^ file in files)
                {
                    songListBox->Items->Add(Path::GetFileName(file));
                }
                
                if (files->Length == 0)
                {
                    statusLabel->Text = String::Format("Status: No songs found in {0}", MUSIC_FOLDER);
                }
                else
                {
                    statusLabel->Text = String::Format("Status: Loaded {0} songs", files->Length);
                }
            }
            catch (Exception^ ex)
            {
                MessageBox::Show("Error loading songs: " + ex->Message, "Error",
                    MessageBoxButtons::OK, MessageBoxIcon::Error);
            }
        }

        void EnsureDirectoriesExist()
        {
            if (!Directory::Exists(MUSIC_FOLDER))
            {
                Directory::CreateDirectory(MUSIC_FOLDER);
            }
            if (!Directory::Exists(MODELS_FOLDER))
            {
                Directory::CreateDirectory(MODELS_FOLDER);
            }
        }

        void StartListening_Click(Object^ sender, EventArgs^ e)
        {
            if (!isListening)
            {
                isListening = true;
                startListeningButton->Text = L"Stop Listening";
                statusLabel->Text = L"Status: Listening for commands...";
                recognitionTimer->Start();
            }
            else
            {
                isListening = false;
                startListeningButton->Text = L"Start Listening";
                statusLabel->Text = L"Status: Ready";
                recognitionTimer->Stop();
            }
        }

        void RecognitionTimer_Tick(Object^ sender, EventArgs^ e)
		{
			recognitionTimer->Stop();
			try
			{
				Console::WriteLine("Starting recording...");
				statusLabel->Text = "Status: Starting recording...";
				recorder->StartRecording();
        
				Console::WriteLine("Processing recording...");
				statusLabel->Text = "Status: Processing recording...";
				recorder->SaveToFile();
				recorder->ProcessRecording();
				recorder->GenerateObservationSequence();
        
				Console::WriteLine("Recognizing command...");
				statusLabel->Text = "Status: Recognizing command...";
				String^ command = RecognizeCommand();
				commandLabel->Text = "Last Command: " + command;
                
                if (command == "create")
                {
                    CreateNewPlaylist();
                }
                else if (command == "add")
                {
                    Sleep(1000);
                    recorder->StartRecording();
                    recorder->SaveToFile();
                    recorder->ProcessRecording();
                    recorder->GenerateObservationSequence();
                    
                    int songIndex = RecognizeDigit();
                    AddSongToPlaylist(songIndex);
                }
                else if (command == "play")
                {
                    PlayCurrentSong();
                }
                else if (command == "stop")
                {
                    StopPlayback();
                }
                else if (command == "shuffle")
                {
                    ShufflePlaylist();
                }
            }
            catch (Exception^ ex)
			{
				Console::WriteLine("Error in RecognitionTimer_Tick: " + ex->Message);
				Console::WriteLine("Stack Trace: " + ex->StackTrace);
				statusLabel->Text = "Status: Error processing command";
			}
    
			if (isListening)
			{
				recognitionTimer->Start();
			}
        }

        void CreateNewPlaylist()
        {
            String^ playlistName = String::Format("Playlist_{0}", playlists->Count + 1);
            playlists->Add(playlistName);
            currentPlaylist = playlistName;
            currentPlaylistLabel->Text = "Current Playlist: " + playlistName;
            playlistListBox->Items->Clear();
            statusLabel->Text = "Status: Created new playlist";
        }

        void AddSongToPlaylist(int index)
        {
            if (currentPlaylist == nullptr)
            {
                statusLabel->Text = "Status: No playlist selected";
                return;
            }
            
            if (index >= 0 && index < songListBox->Items->Count)
            {
                String^ song = safe_cast<String^>(songListBox->Items[index]);
                playlistListBox->Items->Add(song);
                statusLabel->Text = String::Format("Status: Added song {0}", song);
            }
            else
            {
                statusLabel->Text = "Status: Invalid song index";
            }
        }

        void PlayCurrentSong()
        {
            if (playlistListBox->SelectedItem != nullptr)
            {
                String^ songPath = Path::Combine(MUSIC_FOLDER,
                    safe_cast<String^>(playlistListBox->SelectedItem));
                
                if (File::Exists(songPath))
                {
                    try
                    {
                        mediaPlayer->SoundLocation = songPath;
                        mediaPlayer->Play();
                        nowPlayingLabel->Text = "Now Playing: " + 
                            safe_cast<String^>(playlistListBox->SelectedItem);
                        statusLabel->Text = "Status: Playing " + 
                            safe_cast<String^>(playlistListBox->SelectedItem);
                    }
                    catch (Exception^ ex)
                    {
                        statusLabel->Text = "Status: Error playing file";
                        MessageBox::Show("Error playing audio: " + ex->Message, "Error",
                            MessageBoxButtons::OK, MessageBoxIcon::Error);
                    }
                }
                else
                {
                    statusLabel->Text = "Status: Song file not found";
                    MessageBox::Show("The selected song file could not be found.", "Error",
                        MessageBoxButtons::OK, MessageBoxIcon::Error);
                }
            }
        }

        void StopPlayback()
        {
            mediaPlayer->Stop();
            nowPlayingLabel->Text = "Now Playing: None";
            statusLabel->Text = "Status: Playback stopped";
        }

        void ShufflePlaylist()
        {
            Random^ rnd = gcnew Random();
            int count = playlistListBox->Items->Count;
            
            while (count > 1)
            {
                count--;
                int k = rnd->Next(count + 1);
                Object^ temp = playlistListBox->Items[k];
                playlistListBox->Items[k] = playlistListBox->Items[count];
                playlistListBox->Items[count] = temp;
            }
            
            statusLabel->Text = "Status: Playlist shuffled";
        }

		void TestButton_Click(System::Object^ sender, System::EventArgs^ e)
		{
			RunTests();
		}
};
}