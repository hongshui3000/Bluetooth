#ifndef CONFIG_H
#define CONFIG_H

/********************Common Parameters**********************/
#define FRAME_SIZE       64
#define FFT_LEVEL        7
#define SAMPLE_RATE      8000
#define DLY_BLOCKS       16
#define BANK_NUM         24

#define PI 3.14159265358979323846


/********************Common Funtion*************************/
#define MAX(a,b) ((a)>(b) ? (a):(b))
#define MIN(a,b) ((a)<(b) ? (a):(b))
#define ABS(x) ((x)>0 ? (x):(-(x)))
#define EXPORT

/********************AEC Parameters*************************/
#define DC_NOTCH_FILTER
#define EMPHASIS_FILTER

#define MDFMethod           1  
#define ScrewedUpNum      50
#define MinLeakage             0.005f
#define ForeUpdateVar1      0.5f
#define ForeUpdateVar2      0.25f
#define BackResetVar           4.0f




/********************NS Parameters*************************/

#define FILTER_BANK

#define GainMethod                    1          /* 0: Loudness domain MMSE  1: log-MMSE*/
#define FilterBankMethod           1          /* 0: Bark domain  1:Mel domian*/
#define SPEECH_PROB_START_DEFAULT     0.35f
#define SPEECH_PROB_CONTINUE_DEFAULT  0.20f
#define NOISE_SUPPRESS_UPDATE_DEFAULT 2.5f
#define NOISE_SUPPRESS_DEFAULT        -20.0f
#define ECHO_SUPPRESS_DEFAULT         -40.0f
#define ECHO_SUPPRESS_ACTIVE_DEFAULT  -15.0f


#endif