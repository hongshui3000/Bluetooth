
#ifndef ECHOCANCEL_H
#define ECHOCANCEL_H

#include "config.h"
#include "CommFun.h"


typedef struct 
{

	/* basic info*/
	int iDlyBlocks;
	int iFrameSize;
	int iFFTLevel;
	int iSampleRate;


#ifdef DC_NOTCH_FILTER
		float fDCNotchRadius;
	    float *fDCNotchMem; 
#endif


#ifdef EMPHASIS_FILTER
		float fEmphasisCoe;
	    float *fEmphasisMem; 
#endif

	int iFrameNum;
	int iAdaptedCon;
	float fAdaptedSum;
	int iSaturated;
	int iScrewedUp;

	float fLeakEstimate;

	float Pey;
	float Pyy;
	float Davg1;
	float Davg2;
	float Dvar1;
	float Dvar2;
  
	float *d;     
	float *x; 
	float *X;
	float *Y;	
	float *E;
	float *Fore;
	float *Back;
	float *BackProp;
	float *BackTemp;
    float *e;
    float *y;
	float *y_overlap;
    float *PHI; 

	float *Xf;
	float *Yf;
	float *Ef;
	float *Yh;
	float *Eh;
	
	float *power;
	float *power_1;
		
	float *window;

}EchoCancelState;

EchoCancelState *EchoCancelStateInit(int iDlyBlocks, int iFrameSize,int  iFFTLevel, int iSampleRate);

void EchoCancelStateDestory(EchoCancelState *st);

void EchoCancelRun(EchoCancelState *st, float *far_end, float *in, float *far_end_echo, float *echo_residual, float *out);


#endif