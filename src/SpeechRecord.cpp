//Created on Wed Jun 05 16:12:36 2013
//@author: Bill

#include "SpeechRecord.h"

// Set wave format when sampling the audio
int SetWaveFormat(WAVEFORMATEX* wf,int wFormatTag,int nChannels,int nSamplesPerSec,
	int nBlockAlign, int wBitsPerSample, int cbSize)
{
	//int res;

	wf->wFormatTag = wFormatTag;
	wf->nChannels = nChannels;
	wf->nSamplesPerSec = nSamplesPerSec;
	wf->nBlockAlign = nBlockAlign;
	wf->wBitsPerSample = wBitsPerSample;
	wf->cbSize = cbSize;
	wf->nAvgBytesPerSec = nSamplesPerSec * wBitsPerSample / 8;
	
	//IAudioClient *pIAudioClient = NULL;
	//res = pIAudioClient->IsFormatSupported( AUDCLNT_SHAREMODE_SHARED, (const PWAVEFORMATEX)wf, NULL);
	//if ( res != S_OK )
	//{
	//	_debug_print("Format not Supported!",1);
	//	return -1;
	//}
	//else
	//{
	//	_debug_print("Format Supported!");
	//}
	return 0;
}

// Set wave header
int SetWaveHdr(WAVEHDR* waveHeader,DWORD bufLen){
	waveHeader->dwBufferLength = bufLen;
	waveHeader->dwBytesRecorded = bufLen;
	waveHeader->dwUser = 0;
	waveHeader->dwFlags = 0;
	waveHeader->dwLoops = 0;
	waveHeader->lpData = (char *)GlobalLock(GlobalAlloc(GMEM_MOVEABLE|GMEM_SHARE, bufLen));
	return 0;
}

// Open wave input channel
int OpenWaveIn(HWAVEIN* hWaveIn, WAVEFORMATEX* wf)
{
	int res;
	char lpTemp[256];

	res = waveInGetNumDevs();
	if (! res )
	{
		_debug_print("Access WaveIn channel FAILED!",1);
		return -1;
	}
	else
	{
		_debug_print("Access WaveIn channel SUCCEED!");
	}

	// Open wave input channel
	res = waveInOpen(hWaveIn,WAVE_MAPPER, wf, (DWORD)NULL,0L,CALLBACK_WINDOW); 
	if ( res != MMSYSERR_NOERROR )
	{
		sprintf(lpTemp, "Open wave input channel FAILED£¬Error_Code = 0x%x", res );
	    _debug_print(lpTemp,1);
	   return -1;
	}
	else
	{
		_debug_print("Open wave input channel SUCCEED!");
	}
	return 0;
}

// Prepare Wave In Header and allocate memory
int PrepareWaveIn(HWAVEIN* hWaveIn, WAVEHDR* waveHeader, DWORD dataSize)
{
	int res;
	char lpTemp[256];

	waveHeader->dwBufferLength = dataSize;
	waveHeader->dwBytesRecorded = 0;
	waveHeader->dwUser = 0;
	waveHeader->dwFlags = 0;
	waveHeader->dwLoops = 0;
	waveHeader->lpData = (char *)GlobalLock(GlobalAlloc(GMEM_MOVEABLE|GMEM_SHARE, dataSize));
	memset(waveHeader->lpData, 0, dataSize );

	// Prepare Header
	res = waveInPrepareHeader( *hWaveIn, waveHeader, sizeof(WAVEHDR) ); 
	if ( res != MMSYSERR_NOERROR)
	{
		sprintf(lpTemp, "Cannot prepare wave in header£¬Error_Code = 0x%03X", res );
		_debug_print(lpTemp,1);
		return -1;
	}
	else
	{
		_debug_print("Prepare wave in header SUCCEED!");
	}

	res = waveInAddBuffer( *hWaveIn, waveHeader, sizeof(WAVEHDR) );
	if ( res != MMSYSERR_NOERROR) 
	{
		sprintf(lpTemp, "Cannot add buffer for wave in£¬Error_Code = 0x%03X", res );
		_debug_print(lpTemp,1);
		return -1;
	}
	else
	{
		_debug_print("Add buffer for wave in SUCCEED!");
	}
	return 0;
}

// Start recording speech
int StartRecord(HWAVEIN* hWaveIn)
{
	int res;

	res = waveInStart(*hWaveIn);
	if(res != MMSYSERR_NOERROR)
	{
		_debug_print("Start recording FAILED!",1);
		return -1;
	}
	else
	{
		_debug_print("Start recording...",1);
	}
	return 0;
}

// Stop recording speech
int StopRecord(HWAVEIN* hWaveIn, MMTIME* mmTime)
{
	int res;

	res = waveInGetPosition(*hWaveIn, mmTime, sizeof(MMTIME));
	if(res != MMSYSERR_NOERROR)
	{
		_debug_print("Get Position of wave in FAILED!",1);
		return -1;
	}
	else
	{
		_debug_print("Get Position of wave in SUCCEED!");
	}

	res = waveInStop(*hWaveIn);
	if(res != MMSYSERR_NOERROR)
	{
		_debug_print("Stop recording FAILED!",1);
		return -1;
	}
	else
	{
		_debug_print("Stop recording SUCCEED!");
	}

	res = waveInReset(*hWaveIn);
	if(res != MMSYSERR_NOERROR)
	{
		_debug_print("Reset wave in memory FAILED!",1);
		return -1;
	}
	else
	{
		_debug_print("Reset wave in memory SUCCEED!");
	}

	return 0;
}

// Save recorded speech to file
int SaveRecordtoFile(const char* fileName, WAVEFORMATEX* wf, WAVEHDR* waveHeader, MMTIME* mmTime)
{
	//int res;
	DWORD NumToWrite=0;
	DWORD dwNumber = 0;

	if(!waveHeader->dwBytesRecorded)
		waveHeader->dwBytesRecorded = mmTime->u.cb;
	
	HANDLE FileHandle = CreateFile( CString(fileName), GENERIC_WRITE, 
		FILE_SHARE_READ, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);

	dwNumber = FCC("RIFF");
	WriteFile(FileHandle, &dwNumber, 4, &NumToWrite, NULL);

	dwNumber = waveHeader->dwBytesRecorded + 18 + 20;
	WriteFile(FileHandle, &dwNumber, 4, &NumToWrite, NULL);

	dwNumber = FCC("WAVE");
	WriteFile(FileHandle, &dwNumber, 4, &NumToWrite, NULL);

	dwNumber = FCC("fmt ");
	WriteFile(FileHandle, &dwNumber, 4, &NumToWrite, NULL);

	dwNumber = 18L;
	WriteFile(FileHandle, &dwNumber, 4, &NumToWrite, NULL);

	WriteFile(FileHandle, wf, sizeof(WAVEFORMATEX), &NumToWrite, NULL);

	dwNumber = FCC("data");
	WriteFile(FileHandle, &dwNumber, 4, &NumToWrite, NULL);

	dwNumber = waveHeader->dwBytesRecorded;
	WriteFile(FileHandle, &dwNumber, 4, &NumToWrite, NULL);

	WriteFile(FileHandle, waveHeader->lpData, waveHeader->dwBytesRecorded, &NumToWrite, NULL);
	SetEndOfFile( FileHandle );
	CloseHandle( FileHandle );
	FileHandle = INVALID_HANDLE_VALUE;
	
	_debug_print("SaveRecordtoFile SUCCEED!",1);

	return 0;
}

// Release wave in memory
int ReleaseWaveIn(HWAVEIN* hWaveIn, WAVEHDR* waveHeader)
{
	int res;

	res = waveInUnprepareHeader(*hWaveIn, waveHeader, sizeof(WAVEHDR));
	if ( res != MMSYSERR_NOERROR ) 
	{
		_debug_print("UnPrepare Wave In Header FAILED!",1);
		return -1;
	}
	else
	{
		_debug_print("UnPrepare Wave In Header SUCCEED!");
	}

	res = (int)GlobalFree(GlobalHandle( waveHeader->lpData ));
	if ( res != MMSYSERR_NOERROR )
	{
		_debug_print("Global Free FAILED!",1);
		return -1;
	}
	else
	{
		_debug_print("Global Free SUCCEED!");
	}

	return 0;
}

// Close Wave in channel
int CloseWaveIn(HWAVEIN* hWaveIn)
{
	int res;

	res = waveInClose(*hWaveIn);
	if(res != MMSYSERR_NOERROR)
	{
		_debug_print("Close wave in FAILED!",1);
	}
	else
	{
		_debug_print("Close wave in SUCCEED!");
	}
	return 0;
}

// str2num
DWORD FCC(LPSTR lpStr)
{
	DWORD Number = lpStr[0] + lpStr[1] *0x100 + lpStr[2] *0x10000 + lpStr[3] *0x1000000 ;
	return Number;
}

// Print function in debug mode
void _debug_print(const char* content,int flag)
{
	if(flag)
	{
		#ifdef _CONSOLE_APP_
			printf("%s\r\n",content);
		#else _GUI_APP_
			MessageBox(NULL,CString(content),CString("Hint"),MB_OK);
		#endif
	}
	else
	{
		#ifdef _DEBUG_
			#ifdef _CONSOLE_APP_
				printf("%s\r\n",content);
			#else _WINDOWS_APP_
				MessageBox(NULL,CString(content),CString("Hint"),MB_OK);
			#endif
		#endif
	}
}

//**************** End of File **************************