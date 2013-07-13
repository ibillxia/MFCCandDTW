//Created on Wed Jun 05 16:12:36 2013
//@author: Bill

// Debug switch
// You can use '#define _DEBUG_' to start DEBUG mode in your applications
// or without any define to close DEBUG mode.
//#define _DEBUG_

// Application type
#define _CONSOLE_APP_
//#define _WINDOWS_APP_

#include <atlstr.h>
#include <Mmsystem.h>
#include <Audioclient.h>


// Set wave format when sampling the audio
int SetWaveFormat(WAVEFORMATEX* wf,int wFormatTag,int nChannels,int nSamplesPerSec,
	int nBlockAlign, int wBitsPerSample, int cbSize);

// Set wave header
int SetWaveHdr(WAVEHDR* waveHeader,DWORD bufLen);

// Open wave input channel
int OpenWaveIn(HWAVEIN*, WAVEFORMATEX*);

// Prepare Wave In Header and allocate memory
int PrepareWaveIn(HWAVEIN* hWaveIn, WAVEHDR* waveHeader, DWORD dataSize);

// Start recording speech
int StartRecord(HWAVEIN* hWaveIn);

// Stop recording speech
int StopRecord(HWAVEIN* hWaveIn, MMTIME* mmTime);

// Save recorded speech to file
int SaveRecordtoFile(const char* fileName, WAVEFORMATEX* wf, WAVEHDR* waveHeader, MMTIME* mmTime);

// Release wave in memory
int ReleaseWaveIn(HWAVEIN* hWaveIn, WAVEHDR* waveHeader);

// Close Wave in channel
int CloseWaveIn(HWAVEIN* hWaveIn);

// str2num
DWORD FCC(LPSTR lpStr);

// Print function in debug mode
void _debug_print(const char* content,int flag = 0);
//**************** End of File **************************