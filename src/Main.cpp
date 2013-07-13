//************ Main.cpp *****************
//Created on Wed Jun 05 16:12:36 2013
//@author: Bill
// To use it, you should add the dependency of "Winmm.lib" in the project preferences.
//******************************************************************************************//

#include <stdio.h>
#include <math.h>

#include <Windows.h>

#include "FeatureMFCC.h"
#include "SpeechRecord.h"
#pragma comment(lib,"Winmm.lib")

#define MIN2(a,b) ((a) < (b) ? (a) : (b))
#define MIN3(a,b,c) ((MIN2(a,b) < (c)) ? (MIN2(a,b)) : (c))

// Global variables
int tw = 32;
int ts = 16;
int sampleRatio = 16000;
int sampleWidth = 16;
int nBytesPerSample = sampleWidth/8;
int frameSize = int(sampleRatio*tw/1000.0);
int frameStep =  int(sampleRatio*ts/1000.0);
int fmin=(int)((128.0-tw)/ts+1);
int fmax=(int)((768.0-tw)/ts+1);

double simpleVAD(char* waveData,int frameSize);
double dtw(PTRFEAT pf1, PTRFEAT pf2);

double simpleVAD(char* waveData,int frameSize)
{
	int i;
	double vad=0.0,temp=0.0;
	for(i=0;i<frameSize;i++){
		if(nBytesPerSample==2) temp = (double)*((short*)waveData+i);
		else if(nBytesPerSample ==1) temp = (double)*((bool*)waveData+i);
		else if(nBytesPerSample == 4) temp = (double)*((int*)waveData+i);
		if(temp>0) vad += temp;
		else vad -= temp;
	}
	return vad/frameSize;
}

// Dynamic Time Warping
double dtw(PTRFEAT pf1, PTRFEAT pf2){
	int i,j,k,n,len;
	double dist, tsub, tsum;
	double tempdist[4096], dp[4096];

	/* initialize */
	if(pf1->vl != pf2->vl) { fprintf(stderr,"feature vector size not equal!\n"); return -1; }
	n = pf1->vl;
	len = pf1->vs*pf2->vs;
	for(i=0;i<pf1->vs;i++){
		for(j=0;j<pf2->vs;j++){
			tsum = 0.0;  // compute distance
			for(k=0;k<n;k++) {
				tsub = *(pf1->data+i*n+k) - *(pf2->data+j*n+k);
				tsum += tsub*tsub;
			}
			*(tempdist+i*pf1->vs+j) = tsum;
		}
	}

	// initialize for dynamic programing
	for(i=0;i<len;i++){
		*(dp+i) = 0.0;
	}

	/* DTW */
	*dp = *(tempdist);
	for(i=1;i<pf1->vs;i++)
		*(dp+i*pf1->vs) = *(tempdist+i*pf1->vs) + *(dp+(i-1)*pf1->vs);
	for(i=1;i<pf2->vs;i++)
		*(dp+i) = *(tempdist+i) + *(dp+i-1);
	for(i=1;i<pf1->vs;i++){
		for(j=1;j<pf2->vs;j++){
			*(dp+i*pf1->vs+j) = *(tempdist+i*pf1->vs+j) + 
				MIN3(*(dp+(i-1)*pf1->vs+j),*(dp+(i-1)*pf1->vs+j-1),*(dp+i*pf1->vs+j-1) );
		}
	}
	dist = *(dp+(i-1)*pf1->vs+j-1);

	return dist;
}


int main() 
{
	// Initialize for Log
	freopen("log\\error.log","w",stderr);
	FILE *fpvad = fopen("log\\vad.log","w+");
	FILE *fplbl = fopen("log\\vad.lbl","w+");
	FILE *fpres = fopen("log\\cmd.log","w+");
	char str[512];
	SYSTEMTIME st;
	GetLocalTime(&st);
	sprintf(str,"==================================================\n"
		"START TIME: %d-%02d-%02d %02d:%02d:%02d.%03d \n"
		"--------------------------------------------------\n",
		st.wYear,st.wMonth,st.wDay,st.wHour,st.wMinute,st.wSecond,st.wMilliseconds);
	fprintf(fpvad,"%s",str);
	fprintf(fplbl,"%s",str);
	fprintf(fpres,"%s",str);

	// Initialize for MFCC and DTW
	Parameter param;
	Speech sph;
	Feature feat, featleft, featforward, featright;
	sph.fs = sampleRatio;
	sph.nb = sampleWidth;
	feat.data = (double*)malloc(sizeof(double)*16);
	readFeatFromFile(&featleft, "model\\test_left_009.mfcc");
	readFeatFromFile(&featforward, "model\\test_forward_006.mfcc");
	readFeatFromFile(&featright, "model\\test_right_001.mfcc");

	// Initialize for Speech Recording
	WAVEFORMATEX wf;
	HWAVEIN hWaveInAll;
	WAVEHDR waveHdr, waveHdrAll;
	DWORD dataSize = 1000000L;
	DWORD dataSizeAll = dataSize*10;
	MMTIME mmtAll;
	SetWaveFormat(&wf,1,1,sampleRatio,nBytesPerSample,sampleWidth,0);  // Set wave format when sampling the audio
	OpenWaveIn(&hWaveInAll,&wf);   // Open wave input channel
	PrepareWaveIn(&hWaveInAll, &waveHdrAll, dataSizeAll);    // Prepare Wave In Header and allocate memory
	StartRecord(&hWaveInAll);    // Start recording

	// Start Speech Recording, VAD and Speech Recognition.
	int i,j,k,pi;
	DWORD lidx,ridx;
	char curFrame[2048],fn[64],res[32];
	bool sdStart=false;
	double vad,thres = 300.0;
	double d1,d2,d3;
	i=j=pi=1; k=0;
	memset(curFrame,0,sizeof(curFrame));
	while(waveHdrAll.dwBytesRecorded/nBytesPerSample<dataSizeAll-frameSize){
		Sleep(10);
		if(waveHdrAll.dwBytesRecorded/nBytesPerSample>(i-1)*frameStep+frameSize){
			memcpy(curFrame,waveHdrAll.lpData+(i-1)*nBytesPerSample*frameStep,nBytesPerSample*frameSize);
			vad = simpleVAD(curFrame,frameSize);
			fprintf(fpvad,"%lf ",vad);
			if(!(i%10)) fprintf(fpvad,"\n");
			if(vad<thres){k++;if(k>5000/ts)break;}
			else {k=0;}
			if(sdStart == false && vad>thres){
				sdStart = true;pi=i;
				lidx = (i-1)*frameStep*nBytesPerSample;
				printf("\nspeech start index %d: %d\n%.1lf ",j,lidx,vad);
			}else if(sdStart == true && ((vad<thres && i-pi>fmin)||i-pi>fmax)){
				ridx = ((i-2)*frameStep+frameSize)*nBytesPerSample;
				fprintf(fplbl,"%ld %ld speech\n",lidx,ridx);
				printf("%.1lf\nspeech end index %d: %d\tlength:%d\n",vad,j,ridx,ridx-lidx);
				sprintf(fn,"test\\test_%03d.wav",j);
				printf("%s\n",fn);

				// Save Segmented Speech to file
				SetWaveHdr(&waveHdr,ridx-lidx);
				memcpy(waveHdr.lpData,waveHdrAll.lpData+lidx,ridx-lidx);
				waveHdr.dwBytesRecorded = ridx-lidx;
				SaveRecordtoFile(fn,&wf,&waveHdr,NULL);

				// Keywords Recognition via MFCC+DTW
				sph.splen = (ridx-lidx)/nBytesPerSample;
				sph.spdata = (double*)malloc(sizeof(double)*sph.splen);
				if(nBytesPerSample == 2) 
					for(int tti=0;tti<sph.splen;tti++) sph.spdata[tti] = (double)*((short*)waveHdr.lpData+tti);
				else if(nBytesPerSample == 4) 
					for(int tti=0;tti<sph.splen;tti++) sph.spdata[tti] = (double)*((int*)waveHdr.lpData+tti);
				else if(nBytesPerSample == 1) 
					for(int tti=0;tti<sph.splen;tti++) sph.spdata[tti] = (double)*((bool*)waveHdr.lpData+tti);

				// feature extraction and classification
				mfcc(&sph,&feat,&param);
				free(sph.spdata); sph.spdata = NULL;
				
				d1 = dtw(&feat, &featleft);
				d2 = dtw(&feat, &featforward);
				d3 = dtw(&feat, &featright);
				if(d1<d2 && d1<d3 && d1<18000.0) strcpy(res,"left");
				else if(d2<d1 && d2<d3 && d2<18000.0) strcpy(res,"forward");
				else if(d3<d1 && d3<d2 && d3<18000.0) strcpy(res,"right");
				else strcpy(res,"undefined");
				sprintf(fn,"test\\test_%03d_%s.mfcc", j, res);
				saveFeat2File(&feat,fn); // save feature to file

				GetLocalTime(&st);
				fprintf(fpres,"%02d:%02d:%02d.%03d %s\n",st.wHour,st.wMinute,st.wSecond,st.wMilliseconds,res);

				printf("\n\n====================================================\n");
				printf("%lf %lf %lf\n",d1,d2,d3);
				printf("classification results: %s\n",res);
				printf("====================================================\n\n");

				//if(j>14) break;
				sdStart = false;
				pi=i;
				j++;
			}else{
				printf("%.1lf ",vad);
				if(sdStart == true && vad<thres && i-pi<fmin+1){
					ridx = ((i-2)*frameStep+frameSize)*nBytesPerSample;
					fprintf(fplbl,"%ld %ld noise\n",lidx,ridx);
					sdStart = false;
					pi = i;
					printf("\nspeech start index %d: %d CANCELED!\n",j,lidx);
				}
			}
			i++;
		}
	}
	// free memory for MFCC
	clearup();  
	free(feat.data); feat.data = NULL;

	// stop recording
	printf("\n");
	StopRecord(&hWaveInAll,&mmtAll);
	SaveRecordtoFile("test\\myTestAll.wav",&wf,&waveHdrAll,&mmtAll);
	ReleaseWaveIn(&hWaveInAll, &waveHdrAll);
	CloseWaveIn(&hWaveInAll);

	// close log files
	GetLocalTime(&st);
	sprintf(str,"--------------------------------------------------\n"
		"END TIME: %d-%02d-%02d %02d:%02d:%02d.%03d \n"
		"==================================================\n\n\n",
		st.wYear,st.wMonth,st.wDay,st.wHour,st.wMinute,st.wSecond,st.wMilliseconds);
	fprintf(fpvad,"%s",str);
	fprintf(fplbl,"%s",str);
	fprintf(fpres,"%s",str);
	fclose(fpvad);
	fclose(fplbl);
	fclose(fpres);

	system("pause");
	return 0;
}


/*
int main(int argc,char **argv){
	char* GetData(char *pString, long* len);
	int i;
	char *fin,*fout;
	char *speech;
	long len;
	Parameter param;
	Speech sph;
	Feature feat;

	if(argc==3){
		fin = argv[1];
		fout = argv[2];

		feat.data = (double *)malloc(sizeof(double)*16);

		// Initialize for MFCC feature extraction
		sph.fs = sampleRatio;
		sph.nb = sampleWidth;

		speech =  GetData(fin, &len);
		sph.splen = len/2;
		sph.spdata = (double*)malloc(sizeof(double)*sph.splen);
		for(i=0;i<sph.splen;i++)
			sph.spdata[i] = (double)*((short*)speech+i);

		mfcc(&sph,&feat,&param);

		saveFeat2File(&feat,fout);
		free(feat.data);
		free(sph.spdata);
	}
	
	return 0;
}


//获取声音文件数据的函数，pString参数指向要打开的声音文件
char* GetData(char *pString, long* len) 
{
	HMMIO file1;//定义HMMIO文件句柄；
	file1=mmioOpenA(pString,NULL,MMIO_READ);  //以读写模式打开所给的WAVE文件

	long size;
	mmioSeek(file1, 42 ,SEEK_SET);  // 20+sizeof(PCMWAVEFORMAT)
	mmioRead(file1,(char*)&size,4);  //获取WAVE文件的声音数据的大小
	*len = size;

	char *sph;
	sph = (char*)malloc(sizeof(char)*size);  //根据数据的大小申请缓冲区
	mmioSeek(file1, 46, SEEK_SET); //对文件重新定位  24+sizeof(PCMWAVEFORMAT)
	mmioRead(file1,(char*)sph,size);  //读取声音数据
	mmioClose(file1, MMIO_FHOPEN);  //关闭WAVE文件
	return sph;
}

*/

