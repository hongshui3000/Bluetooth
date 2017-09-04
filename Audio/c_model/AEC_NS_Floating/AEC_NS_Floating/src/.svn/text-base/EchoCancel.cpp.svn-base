#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "FileOpera.h"
#include "config.h"
#include "EchoCancel.h"
#include "CommFun.h"


EXPORT EchoCancelState *EchoCancelStateInit(int iDlyBlocks, int iFrameSize, int iFFTLevel, int iSampleRate)
{

	int iCnt1;        

	int M = iDlyBlocks;
    int N = 2*iFrameSize;
	int L = iFFTLevel;
	
	EchoCancelState *st = (EchoCancelState *)malloc(sizeof(EchoCancelState));

	/* basic info*/
	st->iDlyBlocks = iDlyBlocks;
	st->iFrameSize = iFrameSize;
	st->iFFTLevel = iFFTLevel;
	st->iSampleRate = iSampleRate;


	st->iAdaptedCon = 0;
	st->fAdaptedSum = 0;
	st->iSaturated = 0;
	st->iFrameNum = 0;
	st->iScrewedUp = 0;
	
	st->fLeakEstimate = 0.0;

	st->Pey = 1.0f;
	st->Pyy = 1.0f;
	st->Davg1 = st->Davg2 = st->Dvar1 = st->Dvar2 = 0.0;
	/* DSP-related parameters*/

#ifdef DC_NOTCH_FILTER
	if(iSampleRate < 12e3)
	{
         st->fDCNotchRadius = 0.9f;
	}
	else if(iSampleRate < 24e3)
	{
         st->fDCNotchRadius = 0.982f;
	}
	else
	{
         st->fDCNotchRadius = 0.992f;
	}

	st->fDCNotchMem = (float *)malloc(2*sizeof(float)); 
	st->fDCNotchMem[0] = st->fDCNotchMem[1] = 0.0;
#endif


#ifdef EMPHASIS_FILTER
	st->fEmphasisCoe = 0.9f;
	st->fEmphasisMem = (float *)malloc(3*sizeof(float)); 
	st->fEmphasisMem[0] = st->fEmphasisMem[1] = st->fEmphasisMem[2] = 0.0;
#endif


	st->d = (float *)malloc(N*sizeof(float));    

	st->x = (float *)malloc(N*sizeof(float)); 
	for(iCnt1=0; iCnt1<N; iCnt1++)
	{
		st->x[iCnt1] = 0.0;
	}

	st->X = (float *)malloc((M+1)*N*sizeof(float));
	for(iCnt1=0; iCnt1<(M+1)*N; iCnt1++)
	{
		st->X[iCnt1] = 0.0;
	}

	st->Y = (float *)malloc(N*sizeof(float));  //initialization at every packet

	st->E = (float *)malloc(N*sizeof(float));
	for(iCnt1=0; iCnt1<N; iCnt1++)
	{
		st->E[iCnt1] = 0.0;
	}

	st->Fore = (float *)malloc(M*N*sizeof(float));
	for(iCnt1=0; iCnt1<M*N; iCnt1++)
	{
		st->Fore[iCnt1] = 0.0;
	}

	st->Back = (float *)malloc(M*N*sizeof(float));
	for(iCnt1=0; iCnt1<M*N; iCnt1++)
	{
		st->Back[iCnt1] = 0.0;
	}

	st->BackProp = (float *)malloc(M*sizeof(float));
	{
		// Ratio of 10 between adaption rate of the first and last block
		float sum = 0.0;
		float decay = exp(-2.4f/M);
		st->BackProp[0] = 0.7f;
		sum = st->BackProp[0];
		for(iCnt1=1; iCnt1<M; iCnt1++)
		{
			st->BackProp[iCnt1] = st->BackProp[iCnt1-1]*decay;
			sum = sum+st->BackProp[iCnt1];
		}
		for(iCnt1=M-1; iCnt1>=0; iCnt1--)
		{
			st->BackProp[iCnt1] = 0.8f*st->BackProp[iCnt1]/sum;
		}
	} 

	st->BackTemp = (float *)malloc(N*sizeof(float));
    st->e = (float *)malloc(N*sizeof(float));
    st->y = (float *)malloc(N*sizeof(float));
	st->y_overlap = (float *)malloc(N*sizeof(float));
    st->PHI = (float *)malloc(N*sizeof(float)); 

	st->Xf = (float *)malloc((iFrameSize+1)*sizeof(float));
	st->Yf = (float *)malloc((iFrameSize+1)*sizeof(float));
	st->Ef = (float *)malloc((iFrameSize+1)*sizeof(float));

	st->Yh = (float *)malloc((iFrameSize+1)*sizeof(float));
	for(iCnt1=0; iCnt1<=st->iFrameSize; iCnt1++)
	{
		st->Yh[iCnt1] = 0.0;
	}

	st->Eh = (float *)malloc((iFrameSize+1)*sizeof(float));
	for(iCnt1=0; iCnt1<=st->iFrameSize; iCnt1++)
	{
		st->Eh[iCnt1] = 0.0;
	}
	
	st->power = (float *)malloc((iFrameSize+1)*sizeof(float));
	for(iCnt1=0; iCnt1<=iFrameSize; iCnt1++)
	{
		st->power[iCnt1] = 0.0;
	}
	
	st->power_1 = (float *)malloc((iFrameSize+1)*sizeof(float));
    for(iCnt1=0; iCnt1<=iFrameSize; iCnt1++)
	{
		st->power_1[iCnt1] = 1.0f;
	}
	
	st->window = (float *)malloc(N*sizeof(float));
    for(iCnt1=0; iCnt1<N; iCnt1++)
	{
		st->window[iCnt1] = 0.5f-0.5f*cos(2*PI*iCnt1/N);
	}

    return st;
}

static void EchoCancelStateReset(EchoCancelState *st)
{
	int iCnt1;
	
    int M = st->iDlyBlocks;
	int N = 2*st->iFrameSize;
	
	st->iFrameNum = 0;
	st->iAdaptedCon = 0;
	st->fAdaptedSum = 0;
	st->iSaturated = 0;
	st->iScrewedUp = 0;

	st->Pey = 1.0f;
	st->Pyy = 1.0f;
	st->Davg1 = st->Davg2 = st->Dvar1 = st->Dvar2 = 0.0;
	
#ifdef DC_NOTCH_FILTER
	st->fDCNotchMem[0] = st->fDCNotchMem[1] = 0.0;
#endif

#ifdef EMPHASIS_FILTER
	st->fEmphasisMem[0] = st->fEmphasisMem[1] = st->fEmphasisMem[2] = 0.0;
#endif
	
	for(iCnt1=0;iCnt1<(M+1)*N; iCnt1++)
	{
		st->X[iCnt1] = 0.0;
	}

	for(iCnt1=0; iCnt1<M*N; iCnt1++)
	{
		st->Back[iCnt1] = 0.0;
		st->Fore[iCnt1] = 0.0;
	}

	for(iCnt1=0; iCnt1<N; iCnt1++)
	{
		st->x[iCnt1] = 0.0;
		st->E[iCnt1] = 0.0;
	}
	
	for(iCnt1=0; iCnt1<=st->iFrameSize; iCnt1++)
	{
		st->power[iCnt1] = 0.0;
		st->power_1[iCnt1] = 1.0f;
		st->Eh[iCnt1] = 0.0;
		st->Yh[iCnt1] = 0.0;
	}
	
	for(iCnt1=0; iCnt1<st->iFrameSize; iCnt1++)
	{
		st->y_overlap[iCnt1] = 0.0;
	}
	 

}

EXPORT void EchoCancelStateDestory(EchoCancelState *st)
{

#ifdef DC_NOTCH_FILTER
	free(st->fDCNotchMem); 
#endif

#ifdef EMPHASIS_FILTER
	free(st->fEmphasisMem);
#endif

	free(st->d);     
	free(st->x); 
	free(st->X);
	free(st->Y);	
	free(st->Fore);
	free(st->Back);
	free(st->BackProp);
	free(st->BackTemp);
    free(st->e);
    free(st->y);
	free(st->y_overlap);
    free(st->PHI); 

	free(st->Xf);
	free(st->Yf);
	free(st->Ef);
	free(st->Yh);
	free(st->Eh);
	
	free(st->power);
	free(st->power_1);	
	free(st->window);
 
	free(st);

}



EXPORT void EchoCancelRun(EchoCancelState *st, float *far_end, float *in, float *far_end_echo, float *echo_residual, float *out)
{
	int iCnt1, iCnt2, iCnt3;
	float fTemp1, fTemp2, fTemp3;

	int M = st->iDlyBlocks;
    int N = 2*st->iFrameSize;
	int L = st->iFFTLevel;


	int iForeUpdate = 0;
	int iBackReset = 0;

	float Sey = 0.0;
	float Syy = 0.0;
	float Sxx = 0.0;
	float Sdd = 0.0;
	float Sff = 0.0;
	float Sbb = 0.0;

	float Dbf = 10.0f;

	float Pey = 1.0f;
	float Pyy = 1.0f;
	float Eh, Yh;

	float fAdaptedRate = 0.0;

	st->iFrameNum++;
	
/*****************************************************************************************/
/*******************************Pre-Process**********************************************/
/*****************************************************************************************/

	 /*Microphone Signal: Remove DC & Pre-Emphasis*/

#ifdef DC_NOTCH_FILTER
    float den2 = st->fDCNotchRadius*st->fDCNotchRadius+0.7f*(1-st->fDCNotchRadius)*(1-st->fDCNotchRadius);
    for(iCnt1=0; iCnt1<st->iFrameSize; iCnt1++)
	{
		float vout = st->fDCNotchMem[0]+in[iCnt1];
		st->fDCNotchMem[0] = st->fDCNotchMem[1]+2*(-in[iCnt1]+st->fDCNotchRadius*vout);
		st->fDCNotchMem[1] = in[iCnt1]-den2*vout;
		in[iCnt1] = st->fDCNotchRadius*vout;
	}
#endif
	

#ifdef EMPHASIS_FILTER
	for(iCnt1=0; iCnt1<st->iFrameSize; iCnt1++)
	{
		fTemp1 = in[iCnt1]-st->fEmphasisCoe*st->fEmphasisMem[0];
		st->fEmphasisMem[0] = in[iCnt1];
		st->d[iCnt1] = fTemp1;
	}
#else
	for(iCnt1=0; iCnt1<st->iFrameSize; iCnt1++)
	{
		st->d[iCnt1] = in[iCnt1];
	}	
#endif
	

	/*Far-end Speech: Pre-Emphasis*/
#ifdef EMPHASIS_FILTER
    for(iCnt1=0; iCnt1<st->iFrameSize; iCnt1++)
	{
		st->x[iCnt1] = st->x[st->iFrameSize+iCnt1];
		fTemp1 = far_end[iCnt1]-st->fEmphasisCoe*st->fEmphasisMem[1];
		st->fEmphasisMem[1] = far_end[iCnt1];
		st->x[st->iFrameSize+iCnt1] = fTemp1;
	}	
#else
	for(iCnt1=0; iCnt1<st->iFrameSize; iCnt1++)
	{
		st->x[iCnt1] = st->x[st->iFrameSize+iCnt1];
		st->x[st->iFrameSize+iCnt1] = far_end[iCnt1];
	}	
#endif
	
	
/*****************************************************************************************/
/*******************************Foreground Output ****************************************/
/*****************************************************************************************/	
	
	/*Shift memory*/
	for(iCnt1=M-1; iCnt1>=0; iCnt1--)
	{
		for(iCnt2=0; iCnt2<N; iCnt2++)
		{
			st->X[(iCnt1+1)*N+iCnt2] = st->X[iCnt1*N+iCnt2];
		}
	}

	FFT(st->x, st->X, N, L);


	/*Get foreground out  (st->e: error Foreground/ leak Foreground)*/
	for(iCnt1=0; iCnt1<N; iCnt1++)
	{
		st->Y[iCnt1] = 0.0;
	}
	for(iCnt1=0; iCnt1<M; iCnt1++)
	{
		st->Y[0] += st->X[iCnt1*N]*st->Fore[iCnt1*N];
		for(iCnt2=1; iCnt2<N-1; iCnt2+=2)
		{
			st->Y[iCnt2] += st->X[iCnt1*N+iCnt2]*st->Fore[iCnt1*N+iCnt2]-st->X[iCnt1*N+iCnt2+1]*st->Fore[iCnt1*N+iCnt2+1];
			st->Y[iCnt2+1] += st->X[iCnt1*N+iCnt2+1]*st->Fore[iCnt1*N+iCnt2]+st->X[iCnt1*N+iCnt2]*st->Fore[iCnt1*N+iCnt2+1];
		}
		st->Y[iCnt2] += st->X[iCnt1*N+iCnt2]*st->Fore[iCnt1*N+iCnt2];
	}
	
	IFFT(st->Y, st->e, N, L);
	
	for(iCnt1=0; iCnt1<st->iFrameSize; iCnt1++)
	{
		st->e[iCnt1] = st->d[iCnt1]-st->e[st->iFrameSize+iCnt1];
	}
	
	/*Relative result (Sff: out Foreground variance)*/
	for(iCnt1=0; iCnt1<st->iFrameSize; iCnt1++)
	{
		Sff += st->e[iCnt1]*st->e[iCnt1];
	}
	
/*****************************************************************************************/
/*******************************Background Update/Output**********************************/
/*****************************************************************************************/		

	/*Adjust proportional adaption rate*/
	if(st->iAdaptedCon == 1)
	{
          fTemp1 = 1.0f;
		  fTemp2 = 1.0f;
		  for(iCnt1=0; iCnt1<M; iCnt1++)
		  {
			  fTemp3 = 1.0;
			  for(iCnt2=0; iCnt2<N; iCnt2++)
			  {
				  fTemp3 += st->Back[iCnt1*N+iCnt2]*st->Back[iCnt1*N+iCnt2];
			  }
			  st->BackProp[iCnt1] = sqrt(fTemp3);
			  
			  if(fTemp1<st->BackProp[iCnt1])
			  {
				  fTemp1 = st->BackProp[iCnt1];
			  }
		  }
		  
		  for(iCnt1=0; iCnt1<M; iCnt1++)
		  {
			  st->BackProp[iCnt1] = st->BackProp[iCnt1]+0.1*fTemp1;
			  fTemp2 = fTemp2+st->BackProp[iCnt1];
		  }
		  
		  for(iCnt1=0; iCnt1<M; iCnt1++)
		  {
			  st->BackProp[iCnt1] = 0.99f*st->BackProp[iCnt1]/fTemp2;
		  }
	}
	
	/*Update background*/
	if(st->iSaturated == 0)
	{
		for(iCnt1=0; iCnt1<M; iCnt1++)
		{
			fTemp1 = st->BackProp[iCnt1]*st->power_1[0];
			st->PHI[0] = fTemp1*st->X[(iCnt1+1)*N]*st->E[0];
			for(iCnt2=1,iCnt3=1; iCnt2<N-1; iCnt2+=2,iCnt3++)
			{
				fTemp1 = st->BackProp[iCnt1]*st->power_1[iCnt3];
				st->PHI[iCnt2] = fTemp1*(st->X[(iCnt1+1)*N+iCnt2]*st->E[iCnt2]+st->X[(iCnt1+1)*N+iCnt2+1]*st->E[iCnt2+1]);
				st->PHI[iCnt2+1] = fTemp1*(-st->X[(iCnt1+1)*N+iCnt2+1]*st->E[iCnt2]+st->X[(iCnt1+1)*N+iCnt2]*st->E[iCnt2+1]);
			}
			fTemp1 = st->BackProp[iCnt1]*st->power_1[iCnt3];
			st->PHI[iCnt2] = fTemp1*st->X[(iCnt1+1)*N+iCnt2]*st->E[iCnt2];
			
			for(iCnt2=0; iCnt2<N; iCnt2++)
			{
				st->Back[iCnt1*N+iCnt2] += st->PHI[iCnt2];
			}

		}
	}
	else
	{
		st->iSaturated --;
	}

	// Saturated
	for(iCnt1=0; iCnt1<st->iFrameSize; iCnt1++)
	{
		if(in[iCnt1]<=-32000 || in[iCnt1]>=32000)
		{
			if(st->iSaturated == 0)
			{
				st->iSaturated = 1;
			}
		}
	}	
	
	/*MDF /AUMDF  (0: MDF 1:AUMDF)*/
	if(MDFMethod == 0)
	{
		for(iCnt1=0; iCnt1<M; iCnt1++)
		{
			IFFT(st->Back+iCnt1*N, st->BackTemp, N, L);
			for(iCnt2=st->iFrameSize; iCnt2<N; iCnt2++)
			{
				st->BackTemp[iCnt2] = 0.0;
			}
		    FFT(st->BackTemp, st->Back+iCnt1*N, N, L);
		}
	}
	else
	{
		for(iCnt1=0; iCnt1<M; iCnt1++)
		{
			if(iCnt1==0 || st->iFrameNum%(M-1)==iCnt1-1)
			{
				IFFT(st->Back+iCnt1*N, st->BackTemp, N, L);

			    for(iCnt2=st->iFrameSize; iCnt2<N; iCnt2++)
			    {
				    st->BackTemp[iCnt2] = 0.0;
			    }
			    FFT(st->BackTemp, st->Back+iCnt1*N, N, L);	
			}
		}
	}

    /*Get Background out  (st->y: ~/ leak Background)*/
	for(iCnt1=0; iCnt1<N; iCnt1++)
	{
		st->Y[iCnt1] = 0.0;
	}
	for(iCnt1=0; iCnt1<M; iCnt1++)
	{
		st->Y[0] += st->X[iCnt1*N]*st->Back[iCnt1*N];
		for(iCnt2=1; iCnt2<N-1; iCnt2+=2)
		{
			st->Y[iCnt2] += st->X[iCnt1*N+iCnt2]*st->Back[iCnt1*N+iCnt2]-st->X[iCnt1*N+iCnt2+1]*st->Back[iCnt1*N+iCnt2+1];
			st->Y[iCnt2+1] += st->X[iCnt1*N+iCnt2+1]*st->Back[iCnt1*N+iCnt2]+st->X[iCnt1*N+iCnt2]*st->Back[iCnt1*N+iCnt2+1];
		}
		st->Y[iCnt2] += st->X[iCnt1*N+iCnt2]*st->Back[iCnt1*N+iCnt2];
	}

	IFFT(st->Y, st->y, N, L);

	/*Difference of background and foreground out*/
	for(iCnt1=0; iCnt1<st->iFrameSize; iCnt1++)
	{
		st->e[iCnt1] = st->e[st->iFrameSize+iCnt1]-st->y[st->iFrameSize+iCnt1];
		Dbf += st->e[iCnt1]*st->e[iCnt1];
	}
	
	/*Get Background out  (st->e: error Bckground/ leak Foreground)*/
	for(iCnt1=0; iCnt1<st->iFrameSize; iCnt1++)
	{
		st->e[iCnt1] = st->d[iCnt1]-st->y[st->iFrameSize+iCnt1];
	}
	for(iCnt1=0; iCnt1<st->iFrameSize; iCnt1++)
	{
		Sbb += st->e[iCnt1]*st->e[iCnt1];
	}
	
	
/*****************************************************************************************/
/***********************Foreground update / Background reset *****************************/
/*****************************************************************************************/		
	st->Davg1 = 0.6f*st->Davg1+0.4f*(Sff-Sbb);
	st->Davg2 = 0.85f*st->Davg2+0.15f*(Sff-Sbb);
	st->Dvar1 = 0.36f*st->Dvar1+0.16f*Sff*Dbf;
	st->Dvar2 = 0.7225f*st->Dvar2+0.0225f*Sff*Dbf;
	
	/*Logic for updating the foreground filter*/
	if((Sff-Sbb)*ABS(Sff-Sbb) > Sff*Dbf)
	{
		iForeUpdate = 1;
	}
	else if(st->Davg1*ABS(st->Davg1) > ForeUpdateVar1*st->Dvar1 )
	{
		iForeUpdate = 1;
	}
	else if(st->Davg2*ABS(st->Davg2) > ForeUpdateVar2*st->Dvar2)
	{
		iForeUpdate = 1;
	}


	if(iForeUpdate == 1)
	{
		/*Update foreground*/
		st->Davg1 = st->Davg2 = st->Dvar1 = st->Dvar2 = 0.0;
		for(iCnt1=0; iCnt1<M*N; iCnt1++)
		{
			st->Fore[iCnt1] = st->Back[iCnt1];
		}
		/*Apply a smooth transition so as not introduce blocking artifacts*/
		for(iCnt1=0; iCnt1<st->iFrameSize; iCnt1++)
		{
			st->e[st->iFrameSize+iCnt1] = st->window[st->iFrameSize+iCnt1]*st->e[st->iFrameSize+iCnt1]+st->window[iCnt1]*st->y[st->iFrameSize+iCnt1];
		}
	}
	else
	{
		/*Logic for reseting the background filter*/
		if(-(Sff-Sbb)*ABS(Sff-Sbb) > BackResetVar*Sff*Dbf)
		{
			iBackReset = 1;
		}
		else if((-st->Davg1*ABS(st->Davg1)) > BackResetVar*st->Dvar1)
		{
			iBackReset = 1;
		}
		else if((-st->Davg2*ABS(st->Davg2)) > BackResetVar*st->Dvar2)
		{
			iBackReset = 1;
		}

		/*Rest the background filter*/
        if(iBackReset == 1)
		{
			st->Davg1 = st->Davg2 = st->Dvar1 = st->Dvar2 = 0.0;
			Sbb = Sff;
			/*copy foreground filter to background filter*/
			for(iCnt1=0; iCnt1<M*N; iCnt1++)
			{
				st->Back[iCnt1] = st->Fore[iCnt1];
			}
			for(iCnt1=0 ;iCnt1<st->iFrameSize; iCnt1++)
			{
                st->y[st->iFrameSize+iCnt1] = st->e[st->iFrameSize+iCnt1];
				st->e[iCnt1] = st->d[iCnt1]-st->y[st->iFrameSize+iCnt1];
			}
		}
	}

/*****************************************************************************************/
/******************************Final Output***********************************************/
/*****************************************************************************************/		

	// Do some sanity check
	for(iCnt1=0; iCnt1<st->iFrameSize; iCnt1++)
	{
		Sey += st->e[iCnt1]*st->y[st->iFrameSize+iCnt1];
		Syy += st->y[st->iFrameSize+iCnt1]*st->y[st->iFrameSize+iCnt1];
		Sxx += st->x[st->iFrameSize+iCnt1]*st->x[st->iFrameSize+iCnt1];
		Sdd += st->d[iCnt1]*st->d[iCnt1];
	}
	
    // Sxx same with speex(?)
	Sxx = Sxx*2;

	if(!(Syy>=0 && Sxx>=0 && Sbb>=0) || !(Sff<N*1e9 && Syy<N*1e9 && Sxx<N*1e9))
	{
		st->iScrewedUp += ScrewedUpNum;
		for(iCnt1=0; iCnt1<st->iFrameSize; iCnt1++)
		{
			far_end_echo[iCnt1] = st->d[st->iFrameSize+iCnt1];
            out[iCnt1] = 0.0;
		}
	}
	else 
	{
		if(Sff>Sdd+10000*N)
		{
			st->iScrewedUp++;
		}
		else
		{
			st->iScrewedUp = 0;
		}

#ifdef EMPHASIS_FILTER
        for(iCnt1=0; iCnt1<st->iFrameSize; iCnt1++)
        {
			fTemp1 = st->d[iCnt1]-st->e[st->iFrameSize+iCnt1];
			out[iCnt1] = fTemp1+st->fEmphasisCoe*st->fEmphasisMem[2];
			st->fEmphasisMem[2] = out[iCnt1];

			far_end_echo[iCnt1] = in[iCnt1]-out[iCnt1];
	    }	
#else
        for(iCnt1=0; iCnt1<st->iFrameSize; iCnt1++)
        {
			far_end_echo[iCnt1] = st->e[st->iFrameSize+iCnt1];
		    out[iCnt1] = st->d[iCnt1]-far_end_echo[iCnt1];
	    }
#endif
	}
	
/*****************************************************************************************/
/****************************Learning Rate Calculation************************************/
/*****************************************************************************************/	
    //The echo canceller get out of control and reset
	if(st->iScrewedUp >= ScrewedUpNum)
	{
		EchoCancelStateReset(st);
		return;
	}
    /*****************Part1: Leakage Estimation***********/
	//Convert leak and error to frequency domian
	for(iCnt1=0; iCnt1<st->iFrameSize; iCnt1++)
	{
		st->y[iCnt1] = 0;
		st->e[st->iFrameSize+iCnt1] = st->e[iCnt1];
		st->e[iCnt1] = 0;
	}

	FFT(st->e, st->E, N, L);
	FFT(st->y, st->Y, N, L);

    st->Ef[0] = st->E[0]*st->E[0];
	st->Yf[0] = st->Y[0]*st->Y[0];
	for(iCnt1=1,iCnt2=1; iCnt1<N-1; iCnt1+=2,iCnt2++)
	{
		st->Ef[iCnt2] = st->E[iCnt1]*st->E[iCnt1]+st->E[iCnt1+1]*st->E[iCnt1+1];
		st->Yf[iCnt2] = st->Y[iCnt1]*st->Y[iCnt1]+st->Y[iCnt1+1]*st->Y[iCnt1+1];
	}
	st->Ef[iCnt2] = st->E[iCnt1]*st->E[iCnt1];
	st->Yf[iCnt2] = st->Y[iCnt1]*st->Y[iCnt1];

	//Compute filtered spectra and (cross-)correlation
	float beta = (float)st->iFrameSize/st->iSampleRate;
	for(iCnt1=0; iCnt1<=st->iFrameSize; iCnt1++)
	{
		Eh = st->Ef[iCnt1]-st->Eh[iCnt1];
		Yh = st->Yf[iCnt1]-st->Yh[iCnt1];
		Pey = Pey+Eh*Yh;
		Pyy = Pyy+Yh*Yh;
		st->Eh[iCnt1] = (1-beta)*st->Eh[iCnt1]+beta*st->Ef[iCnt1];
		st->Yh[iCnt1] = (1-beta)*st->Yh[iCnt1]+beta*st->Yf[iCnt1];
	}

	Pyy = sqrt(Pyy);
	Pey = Pey/Pyy;

	// Compute corrleation update rate
	fTemp1 = st->iFrameSize/st->iSampleRate;
	float alpha = MIN(2.0f*fTemp1*Syy/Sbb,0.5f*fTemp1);

	st->Pey = (1-alpha)*st->Pey+alpha*Pey;
	st->Pyy = (1-alpha)*st->Pyy+alpha*Pyy;

	if(st->Pyy < 1.0f)
	{
		st->Pyy = 1.0f;
	}

	//We don't really hope to get better than 33dB attenuation anyway

	if(st->Pey > st->Pyy)
	{
		st->Pey = st->Pyy;
	}

	if(st->Pey < MinLeakage*st->Pyy)
	{
		st->Pey = MinLeakage*st->Pyy;
	}

	//Leakage estimation is the linear regression result
	st->fLeakEstimate = st->Pey/st->Pyy;

	/*****************Part2: Residual to Error ratio***********/

	// Residual to error ratio (Average)
	float RER = (0.0001f*Sxx+3.0f*st->fLeakEstimate*Syy)/Sbb;
	if(RER < Sey*Sey/(1+Sbb*Syy))
	{
		RER = Sey*Sey/(1+Sbb*Syy);
	}
	else if(RER > 0.5f)
	{
		RER = 0.5f;
	}


	//Smooth far end speech energy estimate over time
	st->Xf[0] = st->X[0]*st->X[0];
	for(iCnt1=1,iCnt2=1; iCnt1<N-1; iCnt1+=2,iCnt2++)
	{
		st->Xf[iCnt2] = st->X[iCnt1]*st->X[iCnt1]+st->X[iCnt1+1]*st->X[iCnt1+1];
	}
	st->Xf[iCnt2] = st->X[iCnt1]*st->X[iCnt1];
	
	float gamma = 0.35f/M;
	for(iCnt1=0; iCnt1<=st->iFrameSize; iCnt1++)
	{
		st->power[iCnt1] = (1.0-gamma)*st->power[iCnt1]+gamma*st->Xf[iCnt1]+1.0;
	}

	// We consider the filter has had minimal adaption if the following is ture
    if(!st->iAdaptedCon && st->fAdaptedSum>M && st->fLeakEstimate>0.03f)
	{
		st->iAdaptedCon = 1;
	}

	if(st->iAdaptedCon == 1)
	{
		// Normal learning rate calculation once we're past he minimal adaption phase
		for(iCnt1=0; iCnt1<=st->iFrameSize; iCnt1++)
		{
			float r, e;
			//compute frequency-domian adaption mask
			r = st->fLeakEstimate*st->Yf[iCnt1];
			e = st->Ef[iCnt1]+1.0f;
			if(r>0.5f*e)
			{
				r = 0.5f*e;
			}
			r = 0.7f*r+0.3f*RER*e;
			st->power_1[iCnt1] = r/(e*(st->power[iCnt1]+10));
		}
	}
	else
	{
		//Add a small noise to make sure not to have problems when dividing
		Sbb = MAX(Sbb,100*N);
		
		// Temporary adaption rate if filter is not yet adapted enough
	

		if(Sxx > 1000*N)
		{
			fTemp1 = 0.25f*Sxx;
			if(fTemp1 > 0.25f*Sbb)
			{
				fTemp1 = 0.25f*Sbb;
			}
			fAdaptedRate = fTemp1/Sbb;
		}

		for(iCnt1=0 ;iCnt1<=st->iFrameSize; iCnt1++)
		{
			st->power_1[iCnt1] = fAdaptedRate/(st->power[iCnt1]+10.0f);
		}
		st->fAdaptedSum = st->fAdaptedSum+fAdaptedRate;
	}



/****************************************************************************************/
/****************************Save Residual Echo for NS***********************************/
/****************************************************************************************/

	if(st->iAdaptedCon == 1)
	{
		for(iCnt1=0; iCnt1<st->iFrameSize; iCnt1++)
		{
			st->y_overlap[iCnt1] = st->y_overlap[st->iFrameSize+iCnt1];
			st->y_overlap[st->iFrameSize+iCnt1] = in[iCnt1]-out[iCnt1];
		}	

	    /* Apply hanning window*/
	    for(iCnt1=0; iCnt1<N; iCnt1++)
	    {
             st->y[iCnt1] = st->window[iCnt1]*st->y_overlap[iCnt1];
	    }

	    /* Apply FFT*/
	    FFT(st->y, st->Y, N, L);

	    /* Power Spectrum */
	    echo_residual[0] = st->Y[0]*st->Y[0];
	    for(iCnt1=1,iCnt2=1; iCnt1<N-1; iCnt1+=2,iCnt2++)
	    {
		    echo_residual[iCnt2] = st->Y[iCnt1]*st->Y[iCnt1]+st->Y[iCnt1+1]*st->Y[iCnt1+1];
	    }
	    echo_residual[iCnt2] = st->Y[iCnt1]*st->Y[iCnt1];



	    /* Residual Echo */
	    for(iCnt1=0; iCnt1<=st->iFrameSize; iCnt1++)
	    {
		    if(st->fLeakEstimate > 0.5f)
		    {
			    fTemp1 = 1.0f;
		    }
		    else
		    {
			    fTemp1 = 2*st->fLeakEstimate;
		    }

		    echo_residual[iCnt1] = fTemp1*echo_residual[iCnt1];
	    }

	    if(!(echo_residual[0] >=0 && echo_residual[0]<N*1e9))
	    {
		    for(iCnt1=0; iCnt1<=st->iFrameSize; iCnt1++)
		    {
			    echo_residual[iCnt1] = 0.0;
		    }
	    }
	}
	else
	{
		for(iCnt1=0; iCnt1<N; iCnt1++)
		{
			st->y_overlap[iCnt1] = 0.0;
		}	
		for(iCnt1=0; iCnt1<=st->iFrameSize; iCnt1++)
		{
			echo_residual[iCnt1] = 0.0;
		}
	}

	if(0)
	{
		for(iCnt1=0; iCnt1<=st->iFrameSize; iCnt1++)
		{
			printf("ID %d. echo_residual: %f \n", iCnt1,echo_residual[iCnt1]);
		}
	}

}

