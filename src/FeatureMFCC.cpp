/************************************************************************/
/* File: FeatureMFCC.cpp                                                */
/*  This file realizes the task of MFCC feature extraction for          */
/*    time -series like speech signals.                                 */
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

#include "FeatureMFCC.h"

/* global static variables */
static bool _isParamInit = false;
static bool _isWInit = false; // is window function is Initialized
static double *_w;   // weighting window vector
// variables  for FFT
static bool _isFFTInit = false;
static double *_fftbuf = NULL,*_w1c = NULL, *_w3c = NULL;
static long *_jx0 = NULL;
// variables for Filter Bank
static bool _isFBinit = false;
unsigned short *idx;
// variables for DCT
static bool _isDCTInit = false;
static int _dctnin=0, _dctnout = 0;
static double _dctz = 0.0;
static double **_dctk = NULL;
// variables for Lifter
static bool _isLifter = false;
static double *_rlifter = NULL;

// Initialize feature extraction parameters
void _init_param(PTRPARAM param) {
	param->emphco = (double)0.95;		// pre-emphasis coefficient
	param->frameLen = (double)32.0;	// frame length in ms
	param->frameShift = (double)16.0; // frame shift in ms
	param->win = 1;			// weighting window, 1 -> Hamming window

	param->nFilters = 24;		// number of filters in the filter-bank
	param->freqMin = (double)300.0;		// lower frequency bound
	param->freqMax = (double)3400.0;	// higher frequency bound
	param->fftNpts = 512;		// FFT length, should be integer power of 2
	param->_fftm = 9;			// 2^(_fftm-1) < fftNpts <= 2^_fftm

	param->nCeps = 13;	// number of cepstral coefficients
	param->nLifter = 22; // lifter
	param->isInit = true;
	_isParamInit = true;
	return;
}

// Setting feature extraction parameters
void set_param(PTRPARAM param,char* str) {
	char c,*p,s[64];
	int i;

	p = str;
	_init_param(param);   // initialize the *param* as default parameters.
	while(*p){
		if(*p=='-'){
			c=*(++p);   // get parameter type
			while(*(++p)==' ');
			i=0;
			while((*p) && (*p)!=' ') s[i++]=*(p++);   // get parameter value
			s[i] = '\0';

			switch(c){
			case 'k':	// pre-emphasis coefficient
				sscanf(s,"%lf",&(param->emphco));
				break;
			case 'l':		// frame length in ms
				sscanf(s,"%lf",&(param->frameLen));
				break;
			case 'd':	// frame shift in ms
				sscanf(s,"%lf",&(param->frameShift));
				break;
			case 'w':	// weighting window
				sscanf(s,"%d",&(param->win));
				break;
			case 'n':	// number of filters in the filter-bank
				sscanf(s,"%hu",&(param->nFilters));
				break;
			case 'i':		// lower frequency bound
				sscanf(s,"%lf",&(param->freqMin));
				break;
			case 'u':	// higher frequency bound
				sscanf(s,"%lf",&(param->freqMax));
				break;
			case 'b':	// FFT length, should be integer power of 2
				sscanf(s,"%d",&(param->fftNpts));
				break;
			case 'p':	// number of cepstral coefficients
				sscanf(s,"%hu",&(param->nCeps));
				break;
			case 'r':	// liftering
				sscanf(s,"%hu",&(param->nLifter));
				break;
			default:
				break;
			}
		}
		if(!(*p)) break;
		p++;
	}
	param->isInit = true;
	_isParamInit = true;
	return;
}

// Speech pre-emphasis
void _premphasis(PTRSPH sph, double alpha) {
	long long i;
	for(i=sph->splen-1;i>0;i--){
		sph->spdata[i] = sph->spdata[i] - sph->spdata[i-1]*alpha;
	}
	//sph->spdata[0] = (1-alpha)*sph->spdata[0];
	return;
}

// Framing and windowing (frames as rows)
int _vec2frame(PTRSPH sph,PTRFRAME frame,double tw,double ts) {
	int i;
	int nw,ns;
	nw = (int)(sph->fs*tw/1000.0);		// window length
	ns = (int)(sph->fs*ts/1000.0);		// window shift
	frame->fs = nw;
	frame->nf = (int)((sph->splen-nw)/ns+1);

	frame->data = (double*)malloc(sizeof(double)*frame->nf*frame->fs);
	if(!frame->data) {fprintf(stderr,"memory realloc failed in line %d, in file: %s\n",__LINE__,__FILE__);return -1;}
	for(i=0;i<frame->nf;i++){
		memcpy(frame->data+i*nw,sph->spdata+i*ns,sizeof(double)*nw);
	}
	
	return 0;
}

// weighting window
int _weightwin(PTRFRAME frame,int win) {
	int i,j,n;
	
	n=frame->fs;
	if(_isWInit == false){
		_w = (double*)realloc(_w, sizeof(double)*n);
		if(!_w){fprintf(stderr,"memory realloc failed in line %d, in file: %s\n",__LINE__,__FILE__); return -1;}

		switch(win){
		case 1:  // Hamming window
			for(i=0;i<n;i++){ *(_w+i) = (double) (0.54-0.46*(cos(PI2*(double)i/(n)))); }
			break;
		case 2:  // Hanning window
			for(i=0;i<n;i++){ *(_w+i) = (double)(0.5-0.5*(cos(PI2*(double)i/(n-1)))); 	}
			break;
		case 3:  // Blackman window
			for(i=0;i<n;i++){ *(_w+i) = (double)(0.42-0.5*(cos(PI2*(double(i)/(n-1))))+0.08*(cos(PI4*(double(i)/(n-1))))); }
			break;
		default:
			for(i=0;i<n;i++){ *(_w+i) = (double)1.0; 	}
			break;
		}
		_isWInit = true;
	}
	// multiply window
	for(i=0;i<frame->nf;i++){
		for(j=0;j<frame->fs;j++){ 
			*(frame->data+frame->fs*i+j) = *(frame->data+frame->fs*i+j)*_w[j];
		}
	}
	return 0;
}

// FFT initialize: initialize FFT kernel
int _fft_init(int npts,int *_fftm) {
	int m,i,j,ip,nb,lnb,llnb,n2,n4,n6,n8,n12,n16;
	double angd,ang,c,s;

	if(!npts){
		fprintf(stderr,"cannot initialize FFT with 0 points. In line %d, in file: %s\n",__LINE__,__FILE__);
		*_fftm = npts = 0;
		_isFFTInit = false;
		return -1; 
	}

	if(frexp((double)npts,&m) != (double)0.5){
		fprintf(stderr,"FFT N-points must be integer pow of 2.In line %d, in file: %s\n",__LINE__,__FILE__); 
		return -1;
	}
	m--;
	n2 = npts>>1;
	n4 = npts>>2;
	n6 = npts/6;
	n8 = npts>>3;
	n12 = n6>>1;
	n16 = npts>>4;

	_fftbuf = (double*)realloc(_fftbuf, sizeof(double)*npts);
	if(!_fftbuf) {
		fprintf(stderr,"memory realloc failed in line %d, in file: %s\n",__LINE__,__FILE__);
		free(_fftbuf); _fftbuf = NULL;
		_isFFTInit = false;
		return -1;
	}
	_w1c = (double*)realloc(_w1c, sizeof(double)*n4);
	if(!_w1c) { 
		fprintf(stderr,"memory realloc failed in line %d, in file: %s\n",__LINE__,__FILE__);
		free(_fftbuf); _fftbuf = NULL;
		free(_w1c); _w1c = NULL;
		_isFFTInit = false;
		return -1; 
	}
	_w3c = (double*)realloc(_w3c, sizeof(double)*n4);
	if(!_w3c) {
		fprintf(stderr,"memory realloc failed in line %d, in file: %s\n",__LINE__,__FILE__);
		free(_fftbuf); _fftbuf = NULL;
		free(_w1c); _w1c = NULL;
		free(_w3c); _w3c = NULL;
		_isFFTInit = false;
		return -1; 
	}
	_jx0 = (long*) realloc(_jx0, sizeof(long)*npts/3);
	if(!_jx0){
		fprintf(stderr,"memory realloc failed in line %d, in file: %s\n",__LINE__,__FILE__);
		free(_fftbuf); _fftbuf = NULL;
		free(_w1c); _w1c = NULL;
		free(_w3c); _w3c = NULL;
		free(_jx0); _jx0 = NULL;
		_isFFTInit = false;
		return -1; 
	}

	// calculate the values of every points in vector _w1c
	// where _w1c[i] = cos(PI2 * i / npts);
	// the range of the angle is (0,PI/2)
	ang = PI2 / (double)npts;
	angd = PI2 / (double)npts;
	c = cos(angd);
	s = sin(angd);
	_w1c[1] = c;  // cos(PI2/npts);
	_w1c[n4-1] = s; // cos(PI2*(n4-1)/npts);
	_w1c[n8] = 0.707106781186547;  // cos(PI/4);
	for(i=2;i<=n16;i++){ // use trigonometric formula to calculate the remain values
		_w1c[i] = _w1c[i-1]*c - _w1c[n4-i+1]*s;
		_w1c[n4-i] = _w1c[n4-i+1]*c + _w1c[i-1]*s;
		_w1c[n8+i-1] = _w1c[n8+i-2]*c - _w1c[n8-i+2]*s;
		_w1c[n8-i+1] = _w1c[n8-i+2]*c + _w1c[n8+i-2]*s;
	}

	// initialize vector _w3c
	// cosine value of 3 time the angle in vector _w1c
	for(i=1;i<=n12;i++) 
		_w3c[i] = _w1c[i*3];
	for(i=n12+1;i<=n6;i++)
		_w3c[i] = -_w1c[n2-3*i];
	for(i=n6+1;i<n4;i++)
		_w3c[i] = -_w1c[3*i-n2];

	// initialize vector _jx0
	// the purpose of this vector is unknown temporarily
	_jx0[0] = _jx0[1] = _jx0[2] = 0;
	_jx0[3] = n2;
	_jx0[4] = 3*n4;
	ip = 5;
	nb = 3;
	lnb = 1;
	for(i=1;i<=m-4;i++){
		for(j=0;j<nb;j++)
			_jx0[ip+j] = _jx0[ip-nb+j]>>1;
		ip += nb;
		for(j=0;j<lnb;j++){
			_jx0[ip+j] = _jx0[ip-nb-nb-lnb+j] /4 +n2;
			_jx0[ip+j+lnb] = _jx0[ip+j] + n4;
		}
		ip = ip+lnb+lnb;
		llnb = lnb;
		lnb = nb;
		nb = lnb+llnb+llnb;
	}
	*_fftm = m;  // 2^(_fftm-1) < npts <= 2^_fftm;
	_isFFTInit = true;
	return 0;
}

// Rearranges data in the FFT buffer
int _fft_brx(double* x, int m){
	int n, n1, m1, i, ipair, ibr, j, jbr, jbri, k, ia1, ia2, ia3, nh, b;
	double xt;

	n = 1<< m;
	m1 = m>>1;
	n1 = 1<<m1;
	ia1 = n1>>1;
	ia2 = n / n1;
	ia3 = ia1 + ia2;
	nh = n / 2;

	b = (m - m1 - m1) * n1;
	for (ipair = 0; ipair <= b; ipair += n1) {
		ibr = 0;
		xt = x[ipair+ia1];
		x[ipair+ia1] = x[ipair+ia2];
		x[ipair+ia2] = xt;
		for (i = 1 + ipair; i < ia1 + ipair; i++) {
			k = nh;
			if (k <= ibr)
				do { ibr -= k; k = k/2; } while (k <= ibr);
				ibr += k;
				xt = x[ibr+i+ia1];
				x[ibr+i+ia1] = x[ibr+i+ia2];
				x[ibr+i+ia2] = xt;
				jbr = 0;

				if (m < 4) continue;

				for (j = ibr + ipair; j < ibr + i; j++) {
					jbri = jbr + i;
					xt = x[jbri];
					x[jbri] = x[j];
					x[j] = xt;
					xt = x[jbri+ia1];
					x[jbri+ia1] = x[j+ia2];
					x[j+ia2] = xt;
					xt = x[jbri+ia2];
					x[jbri+ia2] = x[j+ia1];
					x[j+ia1] = xt;
					xt = x[jbri+ia3];
					x[jbri+ia3] = x[j+ia3];
					x[j+ia3] = xt;
					k = nh;
					if(k <= jbr)
						do { jbr -= k; k = k/2; } while (k <= jbr);
						jbr += k;
				}
		}
	}
	return 0;
}

// FFT main body
int _fft_mbd(double* x, int m) {
	int i, i0, i1, i2, i3, i4, i5, i6, i7, ib, istep, ia0, ia1, ia2, ia3;
	int n, ib0, ib1, ib2, ib3, j, jstep, n2, n4 ,n8, nd4, nb, lnb, llnb, k, sgn;
	double c2, c3, d2, d3, r1, r2, r3, r4, t0, t1, t2;
	const double rac2s2 = 0.707106781186547;

	n = 1 << m;
	nd4 = n>>2;
	sgn = ((m%2) == 0) ? 1 : -1;
	nb = (n / 2 + sgn) / 3;
	lnb = (n - sgn) / 3;
	ib = n / 6;

	for (i = ib; i < ib + nb; i++) {
		i0 = _jx0[i];
		i1 = i0 + 1;
		i2 = i1 + 1;
		i3 = i2 + 1;
		r1 = x[i0] + x[i1];
		t0 = x[i2] + x[i3];
		x[i3] = x[i3] - x[i2];
		x[i1] = x[i0] - x[i1];
		x[i2] = r1 - t0;
		x[i0] = r1 + t0;
	}
	llnb = lnb;
	lnb = nb;
	nb = (llnb - lnb) / 2;
	ib = ib - nb;

	for (i = ib; i < ib + nb; i++) {
		i0 = _jx0[i];
		i4 = i0 + 4;
		i5 = i0 + 5;
		i6 = i0 + 6;
		i7 = i0 + 7;
		r1 = x[i4] - x[i5];
		r3 = x[i4] + x[i5];
		r2 = x[i7] - x[i6];
		r4 = x[i6] + x[i7];
		t0 = r3 + r4;
		x[i6] = r4 - r3;
		x[i4] = x[i0] - t0;
		x[i0] = x[i0] + t0;

		t1 = (r1 + r2) * rac2s2;
		t2 = (r2 - r1) * rac2s2;
		i3 = i0 + 3;
		x[i5] = t2 - x[i3];
		x[i7] = t2 + x[i3];
		i1 = i0 + 1;
		x[i3] = x[i1] - t1;
		x[i1] = x[i1] + t1;
	}

	istep = n / 16;
	n8 = 1;
	n4 = 2;
	n2 = 4;

	for (k = 4; k <= m; k++) {
		llnb = lnb;
		lnb = nb;
		nb = (llnb - lnb) / 2;
		ib = ib - nb;
		n8 = n4;
		n4 = n2;
		n2 = n2 + n2;

		for (i = ib; i < ib + nb; i++) {
			i0 = _jx0[i];
			i1 = i0 + n4;
			i2 = i1 + n4;
			i3 = i2 + n4;
			t0 = x[i2] + x[i3];
			x[i3] = -x[i2] + x[i3];
			x[i2] = x[i0] - t0;
			x[i0] = x[i0] + t0;

			i0 = i0 + n8;
			i1 = i0 + n4;
			i2 = i1 + n4;
			i3 = i2 + n4;
			t1 = (x[i2] - x[i3]) * rac2s2;
			t2 = (x[i2] + x[i3]) * rac2s2;
			x[i2] = -t2 - x[i1];
			x[i3] = -t2 + x[i1];
			x[i1] = x[i0] - t1;
			x[i0] = x[i0] + t1;
		}

		if (n4 < 4) 
			continue;

		for (i = ib; i < ib + nb; i++) {
			jstep = 0;
			for (j = 1; j <= n8 - 1; j++) {
				jstep = jstep + istep;
				ia0 = _jx0[i] + j;

				ia2 = ia0 + n2;
				ib2 = ia2 + n4 - j - j;
				c2 = x[ia2] * _w1c[jstep] + x[ib2] * _w1c[nd4-jstep];
				d2 = -x[ia2] * _w1c[nd4-jstep] + x[ib2] * _w1c[jstep];
				ia3 = ia2 + n4;
				ib3 = ib2 + n4;
				c3 = x[ia3] * _w3c[jstep] - x[ib3] * _w3c[nd4-jstep];
				d3 = x[ia3] * _w3c[nd4-jstep] + x[ib3] * _w3c[jstep];
				ib1 = ia0 + n4;
				t1 = c2 + c3;
				c3 = c2 - c3;
				x[ib2] = -x[ib1] - c3;
				x[ia3] = x[ib1] - c3;
				t2 = d2 - d3;
				ia1 = ib1 - j - j;
				x[ib1] = x[ia1] + t2;
				x[ia1] = x[ia1] - t2;
				d3 = d2 + d3;
				ib0 = ia1 + n4;
				x[ia2] = -x[ib0] + d3;
				x[ib3] = x[ib0] + d3;
				x[ib0] = x[ia0] - t1;
				x[ia0] = x[ia0] + t1;
			}
		}
		istep = istep / 2;
	}
	return 0;
}

// FFT
int _fft(PTRPARAM param) {
	if(_isFFTInit==false){ 
		fprintf(stderr, "Cannot apply FFT before it is correctly Initialized!\n");
		return -1;
	}
	/* for each frame, do FFT */
	_fft_brx(_fftbuf, param->_fftm);		// Rearranges data in the FFT buffer
	_fft_mbd(_fftbuf, param->_fftm);		// FFT main body

	return 0;
}

// Mel idx setting
int _set_mel_idx(unsigned short n, PTRPARAM pp, int fs) {
	unsigned short i;
	double fmin,fmax;
	double f,min,max,d,z;
	
	fmin = pp->freqMin / fs;
	fmax = pp->freqMax / fs;
	idx = (unsigned short*)malloc((n+2)*sizeof(unsigned short));
	if(!idx){ fprintf(stderr,"memory realloc failed in line %d, in file: %s\n",__LINE__,__FILE__); return -1; }

	if(fmax<=fmin) fmax = 0.5;
	if(fmin<0||fmin>0.5||fmax<0||fmax>0.5)	{ 
		fprintf(stderr,"_set_mel_idx(): invalidfrequency range [%lf,%lf] in line %d, in file: %s\n",
			fmin, fmax, __LINE__,__FILE__); 
		free(idx);
		return -1; 
	}

	*idx = (unsigned short)_round(2*fmin*(pp->fftNpts/2-1));
	*(idx+n+1) = (unsigned short)_round(2*fmax*(pp->fftNpts/2-1));
	
	min = _hz2mel(fmin*fs);
	max = _hz2mel(fmax*fs);
	d = (max-min)/(double)(n+1);
	z = (double)(pp->fftNpts/2-1)*2.0/fs;
	f = min;
	for(i=1;i<=n;i++){
		f += d;
		*(idx+i) = (unsigned short)_round(_mel2hz(f)*z);  // index of the original domain
	}
	_isFBinit = true;
	return 0;
}

// filter-bank: Apply triangular filter bank to the energy and return the log of the energy in each band.
int _filterbank(double *buf,int _fftn, unsigned short nfilt, bool powerspec, bool uselog, double *e){
	int i, j, from, to;
	double a, s, m, re, im;

	for(i=0;i<nfilt;i++){
		s = 0.0;
		// ascending step
		from = *(idx+i);
		to = *(idx+i+1);
		a = 1.0/(double)(to-from+1);
		for(j=from;j<to;j++){
			if(j){
				re = *(buf+j);
				im = *(buf+_fftn-j);
				m = powerspec ? (re*re+im*im) : sqrt(re*re+im*im);
			}else{
				m = powerspec ? (*buf* *buf) : fabs(*buf);
			}
			s += m*(1.0-a*(to-j));
		}
		// descending step
		from = to;
		to = *(idx+i+2);
		a = 1.0/(double)(to-from+1);
		for(j=from;j<=to;j++){
			if(j){
				re = *(buf+j);
				im = *(buf+_fftn-j);
				m = powerspec ? (re*re+im*im) : sqrt(re*re+im*im);
			}else{
				m = powerspec ? (*buf * *buf) : fabs(*buf);
			}
			s += m*(1.0-a*(j-from));
		}
		if(uselog) {
			*(e+i) = (s < ENERGY_FLOOR) ? (double)log(ENERGY_FLOOR) : (double)log(s);
		}else{
			*(e+i) = (s < ENERGY_FLOOR) ? (double)ENERGY_FLOOR : (double)s;
		}
	}
	return 0;
}

// DCT initialize
int _dct_init(unsigned short nin, unsigned short nout){
	double *kp;
	unsigned short i,j;
	if(nin && nout){
		double **p = _dctk;
		_dctk = (double**)realloc(_dctk, nout*sizeof(double*));
		if(!_dctk){
			fprintf(stderr,"memory realloc failed in line %d, in file %s.\n",__LINE__,__FILE__);
			return -1;
		}
		
		for(i=0;i<nout;i++){
			if(!p) _dctk[i] = NULL;
			_dctk[i] = (double*)realloc(_dctk[i],nin*sizeof(double));
			if(!_dctk[i]){
				fprintf(stderr,"memory realloc failed in line %d, in file %s.\n",__LINE__,__FILE__);
				while(i--) free(_dctk[--i]);
				free(_dctk);
				return -1;
			}
			kp = *(_dctk+i);
			for(j=0;j<nin;j++)
				*(kp+j) = (double)cos(PI*(i+1.0)*(j+0.5)/nin);
		}
		_dctz = (double)sqrt(2.0/nin);
		_dctnin = nin;
		_dctnout = nout;
	} else {
		if(_dctk){
			for(i=0;i<_dctnout;i++)
				if(*(_dctk+i)) free(*(_dctk+i));
			free(_dctk);
		}
		_dctnin = _dctnout = 0;
	}
	
	_isDCTInit = true;
	return 0;
}

// DCT
int _dct(double *ip,double *op){
	int i, j;
	double v;
	double *kp;

	if (! _dctnout) {
		fprintf(stderr, "DCT kernel uninitalized");
		return -1;
	}

	for (i = 0; i < _dctnout; i++) {
		kp = *(_dctk+i);
		v = 0.0;
		for (j = 0; j < _dctnin; j++)
			v += ( (*(ip+j)) * (*(kp+j)) );
		*(op+i) = (double)(v * _dctz);
	}
	
	return 0;
}

// liftering: h[i] = 1.0 + l * sin((i + 1) * M_PI / l) / 2.0
int _lifter_set(int l, unsigned short n){
	unsigned short i;
	_rlifter = (double*)malloc(sizeof(double)*n);
	if(!_rlifter){
		fprintf(stderr,"memory alloc failed in line %d, in file %s!\n",__LINE__, __FILE__);
		return -1;
	}
	for(i=0;i<n;i++)
		*(_rlifter+i) = (double)(1.0 + 0.5 * (double)l * sin((double)(i + 1) * PI / (double)l));
	_isLifter = true;
	return 0;
}

/*
** MFCC feature extraction.
** Input: sph - input speech signal
**           param - input feature extraction parameters
** Output: feat - extracted features
**
** Tips: Before calling this function, you should have the *sph* and *param* input initialized and 
**     memory allocated. *param* can be initialized by *set_param* function, if it havn't been initialized, 
**     we will call the *init_param* function to initialize it with a group of default parameters. 
*/
int mfcc(PTRSPH sph, PTRFEAT feat, PTRPARAM param) {
	int i,j;
	Frames frame;
	//unsigned short *idx = NULL;
	double *frm = NULL,*e = NULL,*c = NULL;

	/* some initialization */
	// initialize feature extraction parameters
	if(_isParamInit==false){
		_init_param(param);
	}
	// FFT initialize
	if(_isFFTInit==false ){
		if(_fft_init(param->fftNpts, &(param->_fftm))){ // first initialize
			fprintf(stderr,"initialize FFT kernel failed. In line %d, in file: %s\n",__LINE__,__FILE__);
			_isFFTInit = false;
			return -1;
		}
	}
	// filter_bank initialize
	e = (double*)malloc(param->nFilters*sizeof(double));
	if(!e){ 
		fprintf(stderr,"memory realloc failed in line %d, in file: %s\n",__LINE__,__FILE__); 
		free(_fftbuf); free(_w1c); free(_w3c);  free(_jx0);
		return -1; 
	}
	if(_isFBinit == false){
		if( _set_mel_idx(param->nFilters,param,sph->fs)){ 
			fprintf(stderr,"_set_mel_idx failed in line %d, in file: %s\n",  __LINE__,__FILE__); 
			free(_fftbuf); free(_w1c); free(_w3c);  free(_jx0); free(e);
			return -1;
		}
	}
	// DCT initialize
	if(_isDCTInit==false ){
		if(_dct_init(param->nFilters, param->nCeps)){  // first initialize
			fprintf(stderr,"initialize DCT kernel failed. In line %d, in file: %s\n",__LINE__,__FILE__);
			free(_fftbuf); free(_w1c); free(_w3c);  free(_jx0); free(e); free(idx);
			_isDCTInit = false;
			return -1;
		}
	}
	c = (double*)malloc((param->nCeps+1)*sizeof(double));
	if(!c){ 
		fprintf(stderr,"memory realloc failed in line %d, in file: %s\n",__LINE__,__FILE__);
		free(_fftbuf); free(_w1c); free(_w3c);  free(_jx0); free(e); free(idx); free(_dctk);
		return -1; 
	}
	// lifter initialize
	if(_isLifter==false) {
		if(_lifter_set(param->nLifter, param->nCeps)){
			fprintf(stderr,"liftering initialize failed!\n");
			free(_fftbuf); free(_w1c); free(_w3c);  free(_jx0); free(e); free(idx); free(_dctk);free(c);
			return -1;
		}
	}
	
	/* front-end process */
	// pre-emphasis
	_premphasis(sph, param->emphco);
	// divide the signal into frames
	if(_vec2frame(sph, &frame, param->frameLen, param->frameShift)){ return -1; }
	// weighting window
	if(_weightwin(&frame,param->win)){ return -1; }

	// alloc memory for feature results
	feat->vl = param->nCeps;
	feat->vs = frame.nf;
	feat->data = (double*)realloc(feat->data, sizeof(double)*feat->vl*feat->vs);
	if(!feat->data){
		fprintf(stderr,"memory realloc failed for feat->data in line %d, in file %s.\n",__LINE__,__FILE__);
		free(_fftbuf); free(_w1c); free(_w3c);  free(_jx0); free(e); free(idx); free(_dctk); free(c);
		return -1;
	}

	/* key transformations for each frame */
	frm = frame.data;
	for(i=0;i<frame.nf;i++){ 
		// copy frame to buffer
		memcpy(_fftbuf,frm+i*frame.fs,sizeof(double)*frame.fs);
		for(j=frame.fs;j<param->fftNpts;j++)
			*(_fftbuf+j) = (double)0.0;
		
		if(_fft(param)){  // apply FFT
			fprintf(stderr, "apply FFT failed!\n");
			free(_fftbuf); free(_w1c); free(_w3c);  free(_jx0); free(e); free(idx); free(_dctk); free(c);
			return -1;
		}
		
		if(_filterbank(_fftbuf,param->fftNpts, param->nFilters, 0, 1, e)) {  // apply the filter bank
			fprintf(stderr,"apply filter bank failed\n");
			free(_fftbuf); free(_w1c); free(_w3c);  free(_jx0); free(e); free(idx); free(_dctk); free(c);
			return -1;
		}
		
		if(_dct(e,c)){  // apply DCT
			fprintf(stderr,"apply DCT failed\n");
			free(_fftbuf); free(_w1c); free(_w3c);  free(_jx0); free(e); free(idx); free(_dctk); free(c);
			return -1;
		}

		for(j=0;j<param->nCeps;j++)
			*(c+j) *= *(_rlifter+j);

		/* save the extracted feature to feat->data */
		memcpy(feat->data+feat->vl*i, c, sizeof(double)*feat->vl);
	}

	/* clear up */
	free(e); free(c);free(frame.data);
	return 0;
}

// clear up globally allocated memory
int clearup(void){
	int i;
	free(_w); 
	free(_fftbuf); free(_w1c); free(_w3c); free(_jx0);
	free(idx);
	i=_dctnout;
	while(--i)free(_dctk[i]);
	free(_dctk);
	free(_rlifter);
	return 0;
}

// save Feature Vectors to file
int saveFeat2File(PTRFEAT pf, const char *file) {
	int i,j;
	FILE *fp;
	fp = fopen(file,"w");
	if(!fp){
		fprintf(stderr, "Cannot open file %s\n.",file);
		return -1;
	}
	fprintf(fp,"%d\t%d\n",pf->vs,pf->vl);
	for(i=0;i<pf->vs;i++){
		for(j=0;j<pf->vl;j++)
			fprintf(fp,"%8lf ",*(pf->data+i*pf->vl+j));
		fprintf(fp,"\n");
	}
	fclose(fp);
	return 0;
}

// read feature vectors from file
int readFeatFromFile(PTRFEAT pf, const char *file) {
	int i,j;
	FILE *fp;
	fp = fopen(file,"r");
	fscanf(fp,"%d\t%d",&pf->vs,&pf->vl);
	pf->data = (double*)malloc(sizeof(double)*pf->vs*pf->vl);
	if(!pf->data){
		fprintf(stderr, "Cannot realloc memory for Feature Vectors from file %s.\n",file);
		return -1;
	}
	for(i=0;i<pf->vs;i++){
		for(j=0;j<pf->vl;j++)
			fscanf(fp,"%lf",pf->data+pf->vl*i+j);
	}
	return 0;
}

