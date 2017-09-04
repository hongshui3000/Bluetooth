#ifndef COMMFUN_H
#define COMMFUN_H



typedef struct 
{
   int *bank_left;
   int *bank_right;
   float *filter_left;
   float *filter_right;
   float *scaling;
   int nb_banks;
   int len;
} FilteBankState;


typedef struct
{
	float real;
	float imag;
}compx;

// FilteBankState

FilteBankState *FilteBankInit(int banks, float sampling, int len);

void FilteBankDestory(FilteBankState *bank);

void FilterBankConvert(FilteBankState *bank, float *ps, float *mel);

void FilterBankInverseConvert(FilteBankState *bank, float *mel, float *ps);

// FFT IFFT
void FFT(float *x, float *y, int N, int L);

void IFFT(float *x, float *y, int N, int L);

//Window
float *WelchWindow(int len);

//Hypergeom gain
float HypergeomGain(float xx);




#endif