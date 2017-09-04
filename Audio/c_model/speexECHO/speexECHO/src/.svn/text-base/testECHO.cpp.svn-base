#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "speex_echo.h"
#include "speex_preprocess.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


#define NN 64
#define TAIL 1024


int main()
{
	int iSampleRate = 8000;
	int iSampleNum = 256000;

	FILE *pFarSpeechFid = fopen("farSpeechDec.txt","r");
	FILE *pMicSignalFid = fopen("micSignalDec.txt","r");
	FILE *pNearSpeechFid = fopen("nearSpeechDec.txt","w");


	short *pFarSpeech = (short *)malloc(iSampleNum*sizeof(short));
	short *pMicSignal = (short *)malloc(iSampleNum*sizeof(short));
	short *pNearSpeech = (short *)malloc(iSampleNum*sizeof(short));
	char cTemp1[6];
	char cTemp2[6];
	int iCnt1 = 0;

	for(iCnt1=0; iCnt1<iSampleNum; iCnt1++)
	{
		fscanf(pFarSpeechFid,"%s",&cTemp1);
		fscanf(pMicSignalFid,"%s",&cTemp2);
		pFarSpeech[iCnt1] = atoi(cTemp1);
		pMicSignal[iCnt1] = atoi(cTemp2);
		//printf("%d %d\n",pFarSpeech[iCnt1],pMicSignal[iCnt1]);
	}
    
   // initialization
   SpeexEchoState *st = speex_echo_state_init(NN, TAIL);;
   SpeexPreprocessState *den = speex_preprocess_state_init(NN, iSampleRate);;
   speex_echo_ctl(st, SPEEX_ECHO_SET_SAMPLING_RATE, &iSampleRate);
   speex_preprocess_ctl(den, SPEEX_PREPROCESS_SET_ECHO_STATE, st);

   //remove echo process
   short iFarSpeechReg[NN], iMicSignalReg[NN], iNearSpeechReg[NN];
   int iCnt2 = 0;
   int iCnt3 = 0;
   for (iCnt2=0; iCnt2<iSampleNum; iCnt2=iCnt2+NN)
   {
	   for(iCnt3=0; iCnt3<NN; iCnt3++)
	   {
		   iFarSpeechReg[iCnt3] = pFarSpeech[iCnt2+iCnt3];
		   iMicSignalReg[iCnt3] = pMicSignal[iCnt2+iCnt3];
	   }
       speex_echo_cancellation(st, iMicSignalReg, iFarSpeechReg, iNearSpeechReg);
       speex_preprocess_run(den, iNearSpeechReg);

	   for(iCnt3=0; iCnt3<NN; iCnt3++)
	   {
		   fprintf(pNearSpeechFid,"%d\n",iNearSpeechReg[iCnt3]);
	   }

   }

   speex_echo_state_destroy(st);
   speex_preprocess_state_destroy(den);
   fclose(pFarSpeechFid);
   fclose(pMicSignalFid);
   fclose(pNearSpeechFid);
   return 0;

}