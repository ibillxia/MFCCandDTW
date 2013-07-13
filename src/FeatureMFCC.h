/************************************************************************/
/* File: FeatureMFCC.h                                                  */
/*  This file defines some common macro definitions,                    */
/*  structures, and some function declarations.                         */
/*                                                                      */
/* Reference: SPro 5.0 by Guillaume Gravier (guig@irisa.fr)             */
/*     site: https://gforge.inria.fr/projects/spro/                     */
/*                                                                      */

/*  Copyright (C) 2013 Bill Xia (ibillxia@gmail.com)                    */
/************************************************************************/

/*
* CVS log:
*
* $Author: Bill Xia $
* $Date: 2013/06/25 21:20:12 $
* $Revision: 0.1 $
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <math.h>

#define PI   3.14159265358979323846
#define PI2 6.28318530717958647692
#define PI4 12.56637061435917295384
# define ENERGY_FLOOR 1.0   // floor energy below this threshold

#define _round(x) (((int)ceil(x) - x < x - (int)floor(x)) ? (int)ceil(x) : (int)floor(x))
// Frequency domain to Mel domain
#define _hz2mel(hz) ((double)(1127.0*log(1+((double)hz)/700.0)))
// Mel domain to Frequency domain
#define _mel2hz(mel) ((double)(700.0*(exp((double)(mel)/1127.0)-1.0)))


typedef struct Speech{
	int fs;		// frequency of sampling
	int nb;		// number of bits per sample
	long long splen;	// length of speech
	double* spdata;		// content of speech
}Speech,*PTRSPH;

typedef struct Frames{
	int fs;			// frame size
	int nf;			// number of frames
	double *data;	// frame data, each frame stored as rows
}Frames,*PTRFRAME;

typedef struct Feature{
	int vl;			// feature vector length
	int vs;			// total feature vectors
	double *data;		// feature data, each feature vector stored as rows
}Feature,*PTRFEAT;

typedef struct Parameter{
	bool isInit;			// a flag to identify if the parameters are initialized
	
	double emphco;		// pre-emphasis coefficient
	double frameLen;	// frame length in ms
	double frameShift;	// frame shift in ms
	int win;						// weighting window

	unsigned short nFilters;	// number of filters in the filter-bank
	double freqMin;				// lower frequency bound
	double freqMax;				// higher frequency bound
	int fftNpts;						// FFT length, should be integer power of 2
	int _fftm;							// 

	unsigned short nCeps;	// number of cepstral coefficients  
	unsigned short nLifter; // lifter
}Parameter,*PTRPARAM;

// Initialize feature extraction parameters
void _init_param(PTRPARAM param);
// Setting feature extraction parameters
void set_param(PTRPARAM param,char* str);

// Speech pre-emphasis
void _premphasis(PTRSPH sph, double alpha);
// Framing and windowing (frames as rows)
int _vec2frame(PTRSPH sph,PTRFRAME frame,double tw,double ts);
// weighting window
int _weightwin(PTRFRAME frame,int win);

// FFT initialize: initialize FFT kernel
int _fft_init(int npts,int *_fftm);
// Rearranges data in the FFT buffer
int _fft_brx(double* x, int m);
// FFT main body
int _fft_mbd(double* x, int m);
// FFT
int _fft(PTRPARAM param);

// Mel idx setting
int _set_mel_idx(unsigned short n, PTRPARAM pp, int fs);

// filter-bank: Apply triangular filter bank to the energy and return the log of the energy in each band.
int _filterbank(double *buf,int _fftn, unsigned short nfilt, bool powerspec, bool uselog, double *e);

// DCT initialize
int _dct_init(unsigned short nin, unsigned short nout);
// DCT
int _dct(double *ip,double *op);

// MFCC feature extraction
int mfcc(PTRSPH sph, PTRFEAT feat, PTRPARAM param);

// clear up globally allocated memory
int clearup(void);

// save Feature Vectors to file
int saveFeat2File(PTRFEAT pf, const char *file);

// read feature vectors from file
int readFeatFromFile(PTRFEAT pf, const char *file);

