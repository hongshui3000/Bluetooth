
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "CommFun.h"
#include "FileOpera.h"
#include "NoiseSuppress.h"




EXPORT NoiseSuppressState *NoiseSuppressStateInit(int iFrameSize, int iFFTLevel, int iSampleRate, int iBankNum)
{

	int iCnt1;        

    int N = 2*iFrameSize;
	int L = iFFTLevel;
	int M = iBankNum;
	
	NoiseSuppressState *st = (NoiseSuppressState *)malloc(sizeof(NoiseSuppressState));

	/* basic info*/
	st->iFrameSize = iFrameSize;
	st->iFFTLevel = iFFTLevel;
	st->iSampleRate = iSampleRate;
	st->iBankNum = iBankNum;
	
	st->iFrameNum = 0;                /* Number of frames used for adaption so far*/
	st->iMinNum = 0;                   /* Number of frames used for min search so far*/
	st->iVAD = 0;;                         /* 1:speech 0:noise*/

	/* DSP-related parameters*/
	st->x = (float *)malloc(N*sizeof(float));                /* current frame  time*/
	st->X = (float *)malloc(N*sizeof(float));               /* current frame frequency*/
	st->inBuf = (float *)malloc(iFrameSize*sizeof(float));
	st->outBuf = (float *)malloc(iFrameSize*sizeof(float));
	
	st->window = WelchWindow(N);

	for(iCnt1=0; iCnt1<iFrameSize; iCnt1++)
	{
		st->inBuf[iCnt1] = 0.0;
		st->outBuf[iCnt1] = 0.0;
	}

#ifndef FILTER_BANK

	st->Xf = (float *)malloc((iFrameSize+1)*sizeof(float));             /* current frame power spectrum*/
	st->Nf = (float *)malloc((iFrameSize+1)*sizeof(float));             /* current frame noise power spectrum*/
	st->echoNf = (float *)malloc((iFrameSize+1)*sizeof(float));     /* current frame noise power spectrum*/
	st->Xf_last = (float *)malloc((iFrameSize+1)*sizeof(float));     /* last frame power spectrum*/
	st->priori = (float *)malloc((iFrameSize+1)*sizeof(float));       /* priori SNR*/
	st->porioriAve = (float *)malloc((iFrameSize+1)*sizeof(float));/* Smoothed priori SNR*/
	st->posteriori = (float *)malloc((iFrameSize+1)*sizeof(float)); /* posteriori SNR*/
	st->gainEM = (float *)malloc((iFrameSize+1)*sizeof(float));     /* Ephrain Malah Gain */
	st->gainFloor = (float *)malloc((iFrameSize+1)*sizeof(float));  /* gain floor */
	st->gain = (float *)malloc((iFrameSize+1)*sizeof(float));          /* adjust gain*/
	for(iCnt1=0; iCnt1<=iFrameSize; iCnt1++)
	{
		st->priori[iCnt1] = 0.0;
		st->posteriori[iCnt1] = 1.0f;
		st->porioriAve[iCnt1] = 0.0;
		st->gain[iCnt1] = 1.0f;
		st->Nf[iCnt1] = 0.0;
		st->echoNf[iCnt1] = 0.0;
	}

#else

	st->Xf = (float *)malloc((iFrameSize+1+M)*sizeof(float));                        /* current frame power spectrum*/
	st->Nf = (float *)malloc((iFrameSize+1+M)*sizeof(float));                       /* current frame noise power spectrum*/
	st->echoNf = (float *)malloc((iFrameSize+1+M)*sizeof(float));              /* current frame residual echo noise power spectrum */ 
	st->Xf_last = (float *)malloc((iFrameSize+1+M)*sizeof(float));               /* last frame power spectrum*/
	st->priori = (float *)malloc((iFrameSize+1+M)*sizeof(float));                 /* priori SNR*/
	st->porioriAve = (float *)malloc((iFrameSize+1+M)*sizeof(float));         /* Smoothed priori SNR*/
	st->posteriori = (float *)malloc((iFrameSize+1+M)*sizeof(float));          /* posteriori SNR*/
	st->gainEM = (float *)malloc((iFrameSize+1+M)*sizeof(float));             /* Ephrain Malah Gain */
	st->gainFloor = (float *)malloc((iFrameSize+1+M)*sizeof(float));         /* gain floor*/
	st->gain = (float *)malloc((iFrameSize+1+M)*sizeof(float));                 /* adjust gain*/

	st->bank = FilteBankInit(M,(float)iSampleRate/N,iFrameSize+1);

	for(iCnt1=0; iCnt1<=iFrameSize+M; iCnt1++)
	{
		st->priori[iCnt1] = 0.0;
		st->posteriori[iCnt1] = 1.0f;
		st->porioriAve[iCnt1] = 0.0;
		st->gain[iCnt1] = 1.0f;
		st->Nf[iCnt1] = 0.0;
		st->echoNf[iCnt1] = 0.0;
	}
#endif

	st->s = (float *)malloc((iFrameSize+1)*sizeof(float));                            /* Parameters used for Noise suppression*/
	st->sMin = (float *)malloc((iFrameSize+1)*sizeof(float));                 
	st->sTmp = (float *)malloc((iFrameSize+1)*sizeof(float));
	st->iNoiseUpdate = (int *)malloc((iFrameSize+1)*sizeof(int));             /* probability of speech probability for noise update*/

	for(iCnt1=0; iCnt1<=iFrameSize; iCnt1++)
	{
		st->s[iCnt1] = 0.0;
		st->sMin[iCnt1] = 0.0;
		st->sTmp[iCnt1] = 0.0;
	}

    return st;
}


EXPORT void NoiseSuppressStateDestory(NoiseSuppressState *st)
{

#ifdef FILTER_BANK
	FilteBankDestory(st->bank);
#endif

    free(st->x);              
	free(st->window);  
	free(st->X);              
	free(st->Xf);             

	free(st->inBuf); 
	free(st->outBuf); 

	free(st->Nf);           
	free(st->echoNf);       

	free(st->Xf_last);           

	free(st->priori);          
	free(st->porioriAve);      
	free(st->posteriori);      
	free(st->gainEM);           
	free(st->gainFloor);       
	free(st->gain);           
	free(st->s);             
	free(st->sMin);                 
	free(st->sTmp); 
	free(st->iNoiseUpdate);      
 
	free(st);

}

EXPORT int NoiseSuppressRun(NoiseSuppressState *st, float *in, float *echo_residual, float *out)
{
	int iCnt1, iCnt2;
	float fTemp1, fTemp2, fTemp3;

	int N = 2*st->iFrameSize;
	int L = st->iFFTLevel;
	int M = st->iBankNum;

	st->iFrameNum++;
	st->iMinNum++;

/*****************************************************************************************/
/*****************************Frame Analysis**********************************************/
/*****************************************************************************************/
	/* get FFT input */
	for(iCnt1=0; iCnt1<st->iFrameSize; iCnt1++)
	{
		st->x[iCnt1] = st->inBuf[iCnt1];
		st->x[st->iFrameSize+iCnt1] = in[iCnt1];
		st->inBuf[iCnt1] = in[iCnt1];
	}
    /* apply time domian window*/
	for(iCnt1=0; iCnt1<N; iCnt1++)
	{
		st->x[iCnt1] = st->window[iCnt1]*st->x[iCnt1];
	}
    
	/* FFT */
    FFT(st->x, st->X, N, L);
  
	/* PSD(fisrt FrameSize+1 terms) */
	st->Xf[0] = st->X[0]*st->X[0];
	for(iCnt1=1,iCnt2=1; iCnt1<N-1; iCnt1+=2,iCnt2++)
	{
		st->Xf[iCnt2] = st->X[iCnt1]*st->X[iCnt1]+st->X[iCnt1+1]*st->X[iCnt1+1];
	}
	st->Xf[iCnt2] = st->X[iCnt1]*st->X[iCnt1];

	if(0)
	{
		for(iCnt1=0; iCnt1<=st->iFrameSize; iCnt1++)
		{
			printf("ID %d ->window:%f and frameFreq:%f and PSD:%f\n",iCnt1, st->window[iCnt1], st->X[iCnt1], st->Xf[iCnt1]);
		}
	}
/*****************************************************************************************/
/*****************************Noise Estimate**********************************************/
/*****************************************************************************************/
	/* apply time and frequency domian smoothing*/
	st->s[0] = 0.8f*st->s[0]+0.2f*st->Xf[0];      
	for(iCnt1=1; iCnt1<st->iFrameSize; iCnt1++)
	{
		st->s[iCnt1] = 0.8f*st->s[iCnt1]+0.05f*st->Xf[iCnt1-1]+0.1f*st->Xf[iCnt1]+0.05f*st->Xf[iCnt1+1];
	}
	st->s[iCnt1] = 0.8f*st->s[iCnt1]+0.2f*st->Xf[iCnt1];

	/* the minimum of the past two continuous frames*/
	int minRange;
	if(st->iFrameNum < 100)
	{
		minRange = 15;
	}
	else if(st->iFrameNum < 1000)
	{
		minRange = 50;
	}
	else if(st->iFrameNum < 10000)
	{
		minRange = 150;
	}
	else
	{
		minRange = 300;
	}

	if(st->iMinNum > minRange)
	{
		st->iMinNum = 0.0;
		for(iCnt1=0; iCnt1<=st->iFrameSize; iCnt1++)
		{
			st->sMin[iCnt1] = MIN(st->sTmp[iCnt1], st->s[iCnt1]);
			st->sTmp[iCnt1] = st->s[iCnt1];
		}
	}
	else
	{
		for(iCnt1=0; iCnt1<=st->iFrameSize; iCnt1++)
		{
			st->sMin[iCnt1] = MIN(st->sMin[iCnt1], st->s[iCnt1]);
			st->sTmp[iCnt1] = MIN(st->sTmp[iCnt1], st->s[iCnt1]);
		}
	}

    /* noise update flag */
	for(iCnt1=0; iCnt1<=st->iFrameSize; iCnt1++)
	{
		if(st->s[iCnt1] < st->sMin[iCnt1]*NOISE_SUPPRESS_UPDATE_DEFAULT)
		{
			st->iNoiseUpdate[iCnt1] = 1;
		}
		else
		{
			st->iNoiseUpdate[iCnt1] = 0;
		}
	}


	/* Noise estimate */
	if(st->iFrameNum > 20000)
	{
	     st->iFrameNum = 20000;
	}
    float beta = MAX(0.03f, 1.0f/st->iFrameNum);

	for(iCnt1=0; iCnt1<=st->iFrameSize; iCnt1++)
	{
		if(st->iNoiseUpdate[iCnt1] || st->Xf[iCnt1] < st->Nf[iCnt1])
		{
			st->Nf[iCnt1] = (1-beta)*st->Nf[iCnt1]+beta*st->Xf[iCnt1];
		}
	}

	/* echo noise estimate (from AEC block)*/
	for(iCnt1=0; iCnt1<=st->iFrameSize; iCnt1++)
	{
		st->echoNf[iCnt1] = MAX(0.6f*st->echoNf[iCnt1], echo_residual[iCnt1]);
	}


	if(0)
	{
		for(iCnt1=0; iCnt1<=st->iFrameSize; iCnt1++)
		{
			printf("ID %d ->NoiseUpdate:%d and NOISE: %f and EchoNOISE: %f\n", iCnt1, st->iNoiseUpdate[iCnt1], st->Nf[iCnt1], st->echoNf[iCnt1]);
		}
	}

/*****************************************************************************************/
/*****************************STSA-MMSE***************************************************/
/*****************************************************************************************/

#ifndef FILTER_BANK  /* Gain EM and Speech Prob all in linear Freq*/

/***********************************Priori and Posteriori SNR****************************************/
	if(st->iFrameNum == 1)
	{
		for(iCnt1=0; iCnt1<=st->iFrameSize; iCnt1++)
		{
			st->Xf_last[iCnt1] = st->Xf[iCnt1];
		}
	}

	/* posteriori and poriori SNR estimation*/
	float gamma;
	for(iCnt1=0; iCnt1<=st->iFrameSize; iCnt1++)
	{
		fTemp1 = 1.0f+st->Nf[iCnt1]+st->echoNf[iCnt1];
		st->posteriori[iCnt1] = st->Xf[iCnt1]/fTemp1-1.0f;
		st->posteriori[iCnt1] = MIN(100.0f,st->posteriori[iCnt1]);
		fTemp2 = st->Xf_last[iCnt1]/(st->Xf_last[iCnt1]+fTemp1);
		gamma = 0.1f+0.89f*fTemp2*fTemp2;
		st->priori[iCnt1] = (1.0f-gamma)*st->Xf_last[iCnt1]/fTemp1+gamma*MAX(0.0,st->posteriori[iCnt1]);
		st->priori[iCnt1] = MIN(100.0f,st->priori[iCnt1]);
	}

	if(0)
	{
		for(iCnt1=0; iCnt1<=st->iFrameSize; iCnt1++)
		{
			printf("ID %d ->Posteriori: %f and Priori:%f\n", iCnt1, st->posteriori[iCnt1],st->priori[iCnt1]);
		}
	}

/*****************************************Gain Calculation*******************************************/

	/* apply time and frequency domian smoothing for poriori SNR*/
	st->porioriAve[0] = 0.7f*st->porioriAve[0]+0.3f*st->priori[0];
	for(iCnt1=1; iCnt1<st->iFrameSize; iCnt1++)
	{
		st->porioriAve[iCnt1] = 0.7f*st->porioriAve[iCnt1]+0.075f*st->priori[iCnt1-1]+0.15f*st->priori[iCnt1]+0.075f*st->priori[iCnt1+1];
	}
	st->porioriAve[iCnt1] = 0.7f*st->porioriAve[iCnt1]+0.3f*st->priori[iCnt1];

	/* signal absence parameter*/
	fTemp1 = 0.0f;
	for(iCnt1=0; iCnt1<=st->iFrameSize; iCnt1++)
	{
		fTemp1 = fTemp1+st->porioriAve[iCnt1];
	}
	float PFrame = 0.1f+0.899f*1.0f/(1.0f+0.15f/(fTemp1/(st->iFrameSize+1)));

	/* gain floor from noise and residual echo*/
	float noiseFloor = exp(0.2302585f*NOISE_SUPPRESS_DEFAULT);
	fTemp1 = (1-PFrame)*ECHO_SUPPRESS_DEFAULT+PFrame*ECHO_SUPPRESS_ACTIVE_DEFAULT;
	float echoFloor = exp(0.2302585f*fTemp1);

	/* calculate the gain with signal presence probability*/
	for(iCnt1=0; iCnt1<=st->iFrameSize; iCnt1++)
	{
		float prior_ratio = st->priori[iCnt1]/(1.0f+st->priori[iCnt1]);
		float theta = prior_ratio*(1.0f+st->posteriori[iCnt1]);
		float MM = HypergeomGain(theta);

		st->gainEM[iCnt1] = MIN(1.0f, prior_ratio*MM);
		st->Xf_last[iCnt1] = 0.2f*st->Xf_last[iCnt1]+0.8f*st->gainEM[iCnt1]*st->gainEM[iCnt1]*st->Xf[iCnt1];

		float P1 = 0.199f+0.8f*1.0f/(1.0f+0.15f/(st->porioriAve[iCnt1]));
		float q = 1.0f-PFrame*P1;

		float p = 1.0f/(1.0f+q/(1.0f-q)*(1.0f+st->priori[iCnt1])*exp(-theta));

		st->gainFloor[iCnt1] = sqrt(noiseFloor*st->Nf[iCnt1]+echoFloor*st->echoNf[iCnt1])/sqrt(1.0f+st->Nf[iCnt1]+st->echoNf[iCnt1]);
		st->gainEM[iCnt1] = MAX(st->gainFloor[iCnt1], st->gainEM[iCnt1]);

		if(GainMethod == 0)
		{
			fTemp1 = p*sqrt(st->gainEM[iCnt1])+(1.0f-p)*sqrt(st->gainFloor[iCnt1]);
		    st->gain[iCnt1] = fTemp1*fTemp1;
		}
		else
		{
		    st->gain[iCnt1] = pow(st->gainEM[iCnt1], p)*pow(st->gainFloor[iCnt1], (1.0f-p));
		}
	}


	if(0)
	{
		for(iCnt1=0; iCnt1<=st->iFrameSize; iCnt1++)
		{
			printf("ID %d-> GainEM:%f and Gain:%f\n", iCnt1, st->gainEM[iCnt1], st->gain[iCnt1]);
		}
	}

#else

/***********************************Priori and Posteriori SNR****************************************/
	/* filter bank convert*/
    FilterBankConvert(st->bank, st->Xf, st->Xf+st->iFrameSize+1);
	FilterBankConvert(st->bank, st->Nf, st->Nf+st->iFrameSize+1);
	FilterBankConvert(st->bank, st->echoNf, st->echoNf+st->iFrameSize+1);

	if(0)
	{
		for(iCnt1=0; iCnt1<=st->iFrameSize+M; iCnt1++)
		{
			printf("ID %d ->Xf: %f and Nf:%f and echoNf:%f\n", iCnt1, st->Xf[iCnt1],st->Nf[iCnt1],st->echoNf[iCnt1]);
		}
	}

	if(st->iFrameNum == 1)
	{
		for(iCnt1=0; iCnt1<=st->iFrameSize+M; iCnt1++)
		{
			st->Xf_last[iCnt1] = st->Xf[iCnt1];
		}
	}

	/* posteriori and poriori SNR estimate*/
	float gamma;
	for(iCnt1=0; iCnt1<=st->iFrameSize+M; iCnt1++)
	{
		fTemp1 = 1.0f+st->Nf[iCnt1]+st->echoNf[iCnt1];
		st->posteriori[iCnt1] = st->Xf[iCnt1]/fTemp1-1.0f;
		st->posteriori[iCnt1] = MIN(100.0f,st->posteriori[iCnt1]);
		fTemp2 = st->Xf_last[iCnt1]/(st->Xf_last[iCnt1]+fTemp1);
		gamma = 0.1f+0.89f*fTemp2*fTemp2;
		st->priori[iCnt1] = (1.0f-gamma)*st->Xf_last[iCnt1]/fTemp1+gamma*MAX(0.0,st->posteriori[iCnt1]);
		st->priori[iCnt1] = MIN(100.0f,st->priori[iCnt1]);
	}

	if(0)
	{
		for(iCnt1=0; iCnt1<=st->iFrameSize+M; iCnt1++)
		{
			printf("ID %d ->Posteriori: %f and Priori:%f\n", iCnt1, st->posteriori[iCnt1],st->priori[iCnt1]);
		}
	}


/*****************************************Gain Calculation*******************************************/

	/*apply time domain smoothing of poriori in the filter bank*/
	for(iCnt1=st->iFrameSize+1; iCnt1<=st->iFrameSize+M; iCnt1++)
	{
		st->porioriAve[iCnt1] = 0.7f*st->porioriAve[iCnt1]+0.3f*st->priori[iCnt1];
	}

	/* signal absence parameter*/
	fTemp1 = 0.0;
	for(iCnt1=st->iFrameSize+1; iCnt1<=st->iFrameSize+M; iCnt1++)
	{
		fTemp1 = fTemp1+st->porioriAve[iCnt1];
	}
	float PFrame = 0.1f+0.899f*1.0f/(1.0f+0.15f/(fTemp1/M));

	/* gain floor from noise and residual echo*/
	float noiseFloor = exp(0.2302585f*NOISE_SUPPRESS_DEFAULT);
	fTemp1 = (1-PFrame)*ECHO_SUPPRESS_DEFAULT+PFrame*ECHO_SUPPRESS_ACTIVE_DEFAULT;
	float echoFloor = exp(0.2302585f*fTemp1);

	/* signal presence probability estimate*/
	for(iCnt1=st->iFrameSize+1; iCnt1<=st->iFrameSize+M; iCnt1++)
	{
		float prior_ratio = st->priori[iCnt1]/(1.0f+st->priori[iCnt1]);
		float theta = prior_ratio*(1.0f+st->posteriori[iCnt1]);
		float MM = HypergeomGain(theta);

		st->gainEM[iCnt1] = MIN(1.0f, prior_ratio*MM);

		float P1 = 0.199f+0.8f*1.0f/(1.0f+0.15f/(st->porioriAve[iCnt1]));          //P1 is the poriori probability of signal unabsence
		float q = 1-PFrame*P1;

		st->gain[iCnt1] = 1.0f/(1.0f+q/(1.0f-q)*(1.0f+st->priori[iCnt1])*exp(-theta));   

		st->gainFloor[iCnt1] = sqrt(noiseFloor*st->Nf[iCnt1]+echoFloor*st->echoNf[iCnt1])/sqrt(1+st->Nf[iCnt1]+st->echoNf[iCnt1]);

		st->Xf_last[iCnt1] = 0.2f*st->Xf_last[iCnt1]+0.8f*st->gainEM[iCnt1]*st->gainEM[iCnt1]*st->Xf[iCnt1];
	}

	FilterBankInverseConvert(st->bank, st->gain+st->iFrameSize+1, st->gain);             /* Speech probability */
	FilterBankInverseConvert(st->bank, st->gainEM+st->iFrameSize+1, st->gainEM);         /* Bark Scale Gain*/
	FilterBankInverseConvert(st->bank, st->gainFloor+st->iFrameSize+1, st->gainFloor);    /* Gain floor */

	if(0)
	{
		for(iCnt1=0; iCnt1<=st->iFrameSize+M; iCnt1++)
		{
			printf("ID %d->speech prob:%f and bark gain:%f and gain floor:%f\n",iCnt1,st->gain[iCnt1], st->gainEM[iCnt1], st->gainFloor[iCnt1]);
		}
	}

	/* calculate the gain with signal presence probability */
	for(iCnt1=0; iCnt1<=st->iFrameSize; iCnt1++)
	{
		float prior_ratio = st->priori[iCnt1]/(1+st->priori[iCnt1]);
		float theta = prior_ratio*(1+st->posteriori[iCnt1]);
		float MM = HypergeomGain(theta);
		float p = st->gain[iCnt1];

		float g = MIN(1.0f, prior_ratio*MM);
        st->gainEM[iCnt1] = MIN(g, 3.0f*st->gainEM[iCnt1]); /*Limit with Bark Scale Gain*/

		st->Xf_last[iCnt1] = 0.2f*st->Xf_last[iCnt1]+0.8f*st->gainEM[iCnt1]*st->gainEM[iCnt1]*st->Xf[iCnt1];

		st->gainEM[iCnt1] = MAX(st->gainFloor[iCnt1], st->gainEM[iCnt1]);

		if(GainMethod == 0)
		{
			fTemp1 = p*sqrt(st->gainEM[iCnt1])+(1.0f-p)*sqrt(st->gainFloor[iCnt1]);
		    st->gain[iCnt1] = fTemp1*fTemp1;
		}
		else
		{
		    st->gain[iCnt1] = pow(st->gainEM[iCnt1], p)*pow(st->gainFloor[iCnt1], (1.0f-p));
		}

	}

	if(0)
	{
		for(iCnt1=0; iCnt1<=st->iFrameSize; iCnt1++)
		{
			printf("ID %d->GainEM:%f and Gain:%f\n",iCnt1,st->gainEM[iCnt1], st->gain[iCnt1]);
		}
	}


#endif

/*****************************************************************************************/
/*****************************Frame Analysis**********************************************/
/*****************************************************************************************/
	/* apply gain*/
	st->X[0] = st->gain[0]*st->X[0];
	for(iCnt1=1,iCnt2=1; iCnt1<N-1; iCnt1+=2,iCnt2++)
	{
		st->X[iCnt1] = st->gain[iCnt2]*st->X[iCnt1];
		st->X[iCnt1+1] = st->gain[iCnt2]*st->X[iCnt1+1];
	}
	st->X[iCnt1] = st->gain[iCnt2]*st->X[iCnt1];


	/* apply IFFT*/
	IFFT(st->X, st->x, N, L);

	/* apply time domain window*/
	for(iCnt1=0; iCnt1<N; iCnt1++)
	{
		st->x[iCnt1] = st->x[iCnt1]*st->window[iCnt1];
	}

	for(iCnt1=0; iCnt1<st->iFrameSize; iCnt1++)
	{
		out[iCnt1] = st->x[iCnt1]+st->outBuf[iCnt1];
		st->outBuf[iCnt1] = st->x[st->iFrameSize+iCnt1];
	}

	if(0)
	{
		for(iCnt1=0; iCnt1<st->iFrameSize; iCnt1++)
		{
			printf("ID %d-> in:%f and Out:%f\n", iCnt1, in[iCnt1], out[iCnt1]);
		}
	}

/*****************************************************************************************/
/*****************************VAD Part****************************************************/
/*****************************************************************************************/
	if((PFrame>SPEECH_PROB_START_DEFAULT) || (st->iVAD && PFrame>SPEECH_PROB_CONTINUE_DEFAULT))
	{
		st->iVAD = 1;
	}
	else
	{
		st->iVAD = 0;
	}



	return st->iVAD;

}