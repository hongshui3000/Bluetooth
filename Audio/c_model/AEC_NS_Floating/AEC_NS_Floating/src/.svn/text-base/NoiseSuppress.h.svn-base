
#ifndef NOISESUPPRESS_H
#define NOISESUPPRESS_H

#include "config.h"
#include "CommFun.h"


typedef struct 
{
	/* basic info*/
	int iFrameSize;
	int iFFTLevel;
	int iSampleRate;
	int iBankNum;
	int iFrameNum;                /* Number of frames used for adaption so far*/
	int iMinNum;                  /* Number of frames used for min search so far*/
	int iVAD;                     /* 1:speech 0:noise*/

#ifdef FILTER_BANK 
	FilteBankState *bank;
#endif

	/* DSP-related parameters*/
	float *x;                    /* current frame */
	float *X;                    /* current frame frequency*/
	float *window;               /* time domain window*/
	float *inBuf;                /* input buffer*/
	float *outBuf;               /* output buffer*/
	float *Xf;                   /* current frame PSD*/
	float *Nf;                   /* current noise PSD*/
	float *echoNf;               /* current frame residual echo PSD*/
	float *Xf_last;              /* last frame PSD(signal part)*/
	float *priori;               /* priori SNR*/
	float *porioriAve;           /* smoothed priori SNR*/
	float *posteriori;           /* posteriori SNR*/
	float *gainEM;               /* Ephrain Malah Gain */
	float *gainFloor;            /* gain floor */
	float *gain;                 /* adjust gain*/

	float *s;                    /* smoothed PSD*/
	float *sMin;                 
	float *sTmp;
	int *iNoiseUpdate;           /* noise update*/

}NoiseSuppressState;

NoiseSuppressState *NoiseSuppressStateInit(int iFrameSize,int  iFFTLevel, int iSampleRate, int iBankNum);

int NoiseSuppressRun(NoiseSuppressState *st, float *in, float *echo_residual, float *out);

void NoiseSuppressStateDestory(NoiseSuppressState *st);

#endif