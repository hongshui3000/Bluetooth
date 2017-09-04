#include <stdio.h>
#include <stdlib.h>

#include "FileOpera.h"
#include "config.h"
#include "EchoCancel.h"
#include "NoiseSuppress.h"



void main()
{
	int NN = 256000;

	float *pFarSpeech = (float *)malloc(NN*sizeof(float));
	float *pMicSignal = (float *)malloc(NN*sizeof(float));
	ReadFltFile("farSpeech.txt", NN, pFarSpeech);
	ReadFltFile("micSignal.txt", NN, pMicSignal);
	float *pFarSpeechEcho = (float *)malloc(NN*sizeof(float));
	float *pNearSpeech_AEC = (float *)malloc(NN*sizeof(float));
	float *pNearSpeech_NS = (float *)malloc(NN*sizeof(float));

    // main process
	int iPeriodNum = (int)NN/FRAME_SIZE;
	float *far_end = (float *)malloc(FRAME_SIZE*sizeof(float));
	float *in = (float *)malloc(FRAME_SIZE*sizeof(float));
	float *far_end_echo = (float *)malloc(FRAME_SIZE*sizeof(float));
	float *echo_residual = (float *)malloc((FRAME_SIZE+1)*sizeof(float));
	float *out_aec = (float *)malloc(FRAME_SIZE*sizeof(float));
    float *out_ns = (float *)malloc(FRAME_SIZE*sizeof(float));

	// state initilization
	EchoCancelState *st_aec = EchoCancelStateInit(DLY_BLOCKS, FRAME_SIZE, FFT_LEVEL, SAMPLE_RATE);	
	NoiseSuppressState *st_ns = NoiseSuppressStateInit(FRAME_SIZE, FFT_LEVEL, SAMPLE_RATE, BANK_NUM);	
	int iCnt1;
    int iCnt2;
	for(iCnt1=0; iCnt1<iPeriodNum; iCnt1++)
	{
		for(iCnt2=0; iCnt2<FRAME_SIZE; iCnt2++)
		{
			far_end[iCnt2] = pFarSpeech[iCnt1*FRAME_SIZE+iCnt2];
			in[iCnt2] = pMicSignal[iCnt1*FRAME_SIZE+iCnt2];
		}
		// processing
		EchoCancelRun(st_aec, far_end, in, far_end_echo, echo_residual, out_aec);

		for(iCnt2=0; iCnt2<=FRAME_SIZE; iCnt2++)
		{
			echo_residual[iCnt2] = 0.0;
		}

        NoiseSuppressRun(st_ns, out_aec, echo_residual, out_ns);

		for(iCnt2=0; iCnt2<FRAME_SIZE; iCnt2++)
		{
			pFarSpeechEcho[iCnt1*FRAME_SIZE+iCnt2] = far_end_echo[iCnt2];
		    pNearSpeech_AEC[iCnt1*FRAME_SIZE+iCnt2] = out_aec[iCnt2];
		    pNearSpeech_NS[iCnt1*FRAME_SIZE+iCnt2] = out_ns[iCnt2];
		}
	}

	EchoCancelStateDestory(st_aec);
    NoiseSuppressStateDestory(st_ns);

    WriteFltFile("farSpeechEchoEsti", NN, pFarSpeechEcho);
    WriteFltFile("nearSpeech_aec", NN, pNearSpeech_AEC);
    WriteFltFile("nearSpeech_ns", NN, pNearSpeech_NS);

}