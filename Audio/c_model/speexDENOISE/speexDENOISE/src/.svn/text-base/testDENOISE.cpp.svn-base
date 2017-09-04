#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "speex_preprocess.h"
#include <stdio.h>
#include <stdlib.h>

#define NN 64

int main()
{
	int iSampleRate = 8000;
	int iSampleNum = 256000;

	int iFrameNum = (int)iSampleNum/NN;

	FILE *pNearSpeechFid = fopen("nearSpeechInt.txt","r");
	FILE *pNearSpeechDenoiseFid = fopen("nearSpeechDenoise.txt","w");
	FILE *pVADFid = fopen("VAD.txt","w");

	short *pNearSpeech = (short *)malloc(iSampleNum*sizeof(short));
	char cTemp1[6];
	int iCnt1;
	for(iCnt1=0; iCnt1<iSampleNum; iCnt1++)
	{
		fscanf(pNearSpeechFid,"%s",&cTemp1);
		pNearSpeech[iCnt1] = atoi(cTemp1);
		//printf("%d\n",pNearSpeech[iCnt1]);
	}


// Initialization
   SpeexPreprocessState *st;
   st = speex_preprocess_state_init(NN, iSampleRate);
   int iDENOISEFlg = 1;
   speex_preprocess_ctl(st, SPEEX_PREPROCESS_SET_DENOISE, &iDENOISEFlg);
   int iAGCFlg = 0;
   speex_preprocess_ctl(st, SPEEX_PREPROCESS_SET_AGC, &iAGCFlg);
   int iAGCLevel = 8000;
   speex_preprocess_ctl(st, SPEEX_PREPROCESS_SET_AGC_LEVEL, &iAGCLevel);
   int iDEREVERBFlg = 0;
   speex_preprocess_ctl(st, SPEEX_PREPROCESS_SET_DEREVERB, &iDEREVERBFlg);
   double fDEREVERBDecay = 0.0;
   speex_preprocess_ctl(st, SPEEX_PREPROCESS_SET_DEREVERB_DECAY, &fDEREVERBDecay);
   double fDEREVERBLevel = 0.0;
   speex_preprocess_ctl(st, SPEEX_PREPROCESS_SET_DEREVERB_LEVEL, &fDEREVERBLevel);
   // test
   int iVADFlg = 1;
   speex_preprocess_ctl(st,SPEEX_PREPROCESS_SET_VAD, &iVADFlg);

// Denoise process
   short iNearSpeechTemp[NN];
   int iVAD;
   int iCnt2;
   int iCnt3;
   for(iCnt2=0; iCnt2<iFrameNum; iCnt2++)
   {
	   for(iCnt3=0; iCnt3<NN; iCnt3++)
	   {
		   iNearSpeechTemp[iCnt3] = pNearSpeech[iCnt2*NN+iCnt3];
	   }
	   iVAD = speex_preprocess_run(st, iNearSpeechTemp);
	   for(iCnt3=0; iCnt3<NN; iCnt3++)
	   {
		   fprintf(pNearSpeechDenoiseFid,"%d\n",iNearSpeechTemp[iCnt3]);
	   }
	   fprintf(pVADFid,"%d\n",iVAD);


   }

   speex_preprocess_state_destroy(st);
   return 0;
}
