
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "config.h"
#include "CommFun.h"


FilteBankState *FilteBankInit(int banks, float df, int len)
{
   FilteBankState *bank;
   float max_filterbank, interval_filterbank;
   int i;
   int id1;
   int id2;

   if(FilterBankMethod == 0)
   {
		max_filterbank = 13.1f*atan(.00074f*(df*len))+2.24f*atan((df*len)*(df*len)*1.85e-8f)+1e-4f*(df*len);
   }
   else
   {
	   max_filterbank = 2595.f*log10(1.f+(df*len)/700.f);
   }
   

   interval_filterbank = max_filterbank/(banks-1);
   
   bank = (FilteBankState *)malloc(sizeof(FilteBankState));
   bank->nb_banks = banks;
   bank->len = len;

   bank->bank_left = (int*)malloc(len*sizeof(int));
   bank->bank_right = (int*)malloc(len*sizeof(int));
   bank->filter_left = (float*)malloc(len*sizeof(float));
   bank->filter_right = (float*)malloc(len*sizeof(float));

   /* Think I can safely disable normalisation that for fixed-point (and probably float as well) */

   bank->scaling = (float*)malloc(banks*sizeof(float));

   for (i=0;i<len;i++)
   {
	  float curr_freq;
	  float filterbank;
	  float val;
	  curr_freq = i*df;

      if(FilterBankMethod == 0)
      {
		   filterbank = 13.1f*atan(.00074f*curr_freq)+2.24f*atan(curr_freq*curr_freq*1.85e-8f)+1e-4f*curr_freq;
      }
      else
      {
	       filterbank = 2595.f*log10(1.f+curr_freq/700.f);
      }

      if (filterbank > max_filterbank)
         break;    
      id1 = (int)(floor(filterbank/interval_filterbank));

      if (id1>banks-2)
      {
         id1 = banks-2;
		 val = 1;
      } 
	  else 
	  {
		 val = (filterbank - id1*interval_filterbank)/interval_filterbank;
      }
      id2 = id1+1;
      bank->bank_left[i] = id1;
	  bank->filter_left[i] = 1-val;
      bank->bank_right[i] = id2;
      bank->filter_right[i] = val;

   }
   
   /* Think I can safely disable normalisation for fixed-point (and probably float as well) */
   for (i=0;i<bank->nb_banks;i++)
   {
        bank->scaling[i] = 0;
   }
      
   int id;
   for (i=0;i<bank->len;i++)
   {
      id = bank->bank_left[i];
      bank->scaling[id] += bank->filter_left[i];
      id = bank->bank_right[i];
      bank->scaling[id] += bank->filter_right[i];
   }
   for (i=0;i<bank->nb_banks;i++)
   {
        bank->scaling[i] = 1.0f/bank->scaling[i];
   }
      

   return bank;
} 


void FilteBankDestory(FilteBankState *bank)
{

   free(bank->bank_left);  
   free(bank->bank_right);
   free(bank->filter_left);
   free(bank->filter_right);
   free(bank->scaling);
   free(bank);
}

void FilterBankConvert(FilteBankState *bank, float *ps, float *filterbank)
{
   int i;
   for (i=0;i<bank->nb_banks;i++)
   {
        filterbank[i] = 0;
   }
      

   int id; 
   for (i=0;i<bank->len;i++)
   {
      id = bank->bank_left[i];
	  filterbank[id] += bank->filter_left[i]*ps[i];
      id = bank->bank_right[i];
	  filterbank[id] += bank->filter_right[i]*ps[i];
   }
}

void FilterBankInverseConvert(FilteBankState *bank, float *filterbank, float *ps)
{
   int i;
   float tmp;
   int id1, id2;
   for (i=0;i<bank->len;i++)
   {
      id1 = bank->bank_left[i];
      id2 = bank->bank_right[i];
	  tmp = filterbank[id1]*bank->filter_left[i];
	  tmp += filterbank[id2]*bank->filter_right[i];
	  ps[i] = tmp;
   }
}

/* FFT IFFT */
static compx compxMult(compx a, compx b)
{
	compx c;
	c.real = a.real*b.real-a.imag*b.imag;
	c.imag = a.real*b.imag+a.imag*b.real;
	return c;
};

static compx compxAdd(compx a, compx b)
{
	compx c;
	c.real = a.real+b.real;
	c.imag = a.imag+b.imag;
	return c;
};

static compx compxSub(compx a, compx b)
{
	compx c;
	c.real = a.real-b.real;
	c.imag = a.imag-b.imag;
	return c;
};

static compx compxConj(compx a)
{
	compx c;
	c.real = a.real;
	c.imag = -a.imag;
	return c;
}

static void RaderRank(compx *pData, int N)
{

	compx cTemp;
	int iForwardIdx = 0;
	int iBackwardIdx = 0;
	int iCompar;
	for(iForwardIdx=0; iForwardIdx<N-1; iForwardIdx++)
	{
		if( iForwardIdx<iBackwardIdx)
		{
			cTemp.real = pData[iForwardIdx].real;
			cTemp.imag = pData[iForwardIdx].imag;
			pData[iForwardIdx].real = pData[iBackwardIdx].real;
			pData[iForwardIdx].imag = pData[iBackwardIdx].imag;
			pData[iBackwardIdx].real = cTemp.real;
			pData[iBackwardIdx].imag = cTemp.imag;
		}

		iCompar = N>>1;
		while(iBackwardIdx >= iCompar)
		{
			iBackwardIdx = iBackwardIdx-iCompar;
			iCompar = iCompar>>1;
		}
		iBackwardIdx = iBackwardIdx+iCompar;

	}
}

void FFT(float *x, float *y, int N, int L)
{
	compx *xTemp = (compx *)malloc(N*sizeof(compx));

	int iCnt1, iCnt2, iCnt3;
    for(iCnt1=0; iCnt1<N; iCnt1++)
	{
		xTemp[iCnt1].real = x[iCnt1];
		xTemp[iCnt1].imag = 0;
	}

	RaderRank(xTemp, N);

    int LE, LE1;
	compx u, w, t;

	for(iCnt1=1; iCnt1<=L; iCnt1++)  
	{
        LE = 1<<iCnt1;       
        LE1 = 1<<(iCnt1-1);  
		u.real = 1;
		u.imag = 0;
		w.real = cos(2.0*PI/LE);
		w.imag = -sin(2.0*PI/LE);

		for(iCnt2=1; iCnt2<=LE1; iCnt2++) 
		{
			for(iCnt3=iCnt2-1; iCnt3<=N-1; iCnt3+=LE) 
			{
				t = compxMult(xTemp[iCnt3+LE1], u);
				xTemp[iCnt3+LE1] = compxSub(xTemp[iCnt3],t);
				xTemp[iCnt3] = compxAdd(xTemp[iCnt3],t);
 			}
            u = compxMult(u,w);
		}

	}

	for(iCnt1=0; iCnt1<=N/2; iCnt1++)
	{
		if(iCnt1==0)
		{
			y[iCnt1] = xTemp[iCnt1].real/N;
		}
		else if(iCnt1<N/2)
		{
			y[2*iCnt1-1] = xTemp[iCnt1].real/N;
			y[2*iCnt1] = xTemp[iCnt1].imag/N;
		}
		else
		{
			y[2*iCnt1-1] = xTemp[iCnt1].real/N;
		}

	}

	free(xTemp);

} 

void IFFT(float *x, float *y, int N, int L)
{
	compx *xTemp = (compx *)malloc(N*sizeof(compx));

	int iCnt1, iCnt2, iCnt3;

	for(iCnt1=0; iCnt1<=N/2; iCnt1++)
	{
		if(iCnt1 == 0)
		{
            xTemp[iCnt1].real = x[iCnt1];
			xTemp[iCnt1].imag = 0;
		}
		else if(iCnt1<N/2)
		{
			xTemp[iCnt1].real = x[2*iCnt1-1];
			xTemp[iCnt1].imag = x[2*iCnt1];
			xTemp[N-iCnt1].real = xTemp[iCnt1].real;
			xTemp[N-iCnt1].imag = -xTemp[iCnt1].imag;
		} 
		else
		{
            xTemp[iCnt1].real = x[iCnt1];
			xTemp[iCnt1].imag = 0;
		}		
	}


	RaderRank(xTemp, N);


	int LE, LE1;
	compx u, w, t;

	for(iCnt1=1; iCnt1<=L; iCnt1++)  
	{
        LE = 1<<iCnt1;     
        LE1 = 1<<(iCnt1-1); 
		u.real = 1.0;
		u.imag = 0;
		w.real = cos(2.0*PI/LE);
		w.imag = sin(2.0*PI/LE);

		for(iCnt2=1; iCnt2<=LE1; iCnt2++) 
		{
			for(iCnt3=iCnt2-1; iCnt3<=N-1; iCnt3+=LE) 
			{
				t = compxMult(xTemp[iCnt3+LE1], u);
				xTemp[iCnt3+LE1] = compxSub(xTemp[iCnt3],t);
				xTemp[iCnt3] = compxAdd(xTemp[iCnt3],t);
				xTemp[iCnt3+LE1].real = xTemp[iCnt3+LE1].real/2;
                xTemp[iCnt3+LE1].imag = xTemp[iCnt3+LE1].imag/2;
				xTemp[iCnt3].real = xTemp[iCnt3].real/2;
				xTemp[iCnt3].imag = xTemp[iCnt3].imag/2;
 			}
            u = compxMult(u,w);
		}
	}

	for(iCnt1=0; iCnt1<N; iCnt1++)
	{
		y[iCnt1] = xTemp[iCnt1].real*N;
	}

	free(xTemp);
} 


// Window
 float *WelchWindow(int len)
{
   float *w = (float *)malloc(len*sizeof(float));
   int i;
   for (i=0;i<len;i++)
   {
#if 1                   // enable(1) with hanning window
      float tmp;
	  float x = 4.0f*i/len;
      int inv=0;
      if (x<2.0)
      {
		  x=2.0-x;
         inv=1;
      } 
	  else if (x<3.0)
      {
		 x=x-2.0;
         inv=1;
      } 
	  else
	  {
		  x=4.0-x; /* 4 - x */
      }
	  x = 1.271903*x;
	  tmp = pow(0.5-0.5*cos(0.5*PI*x),2);
      if (inv)
		 tmp = 1-tmp;
	  w[i]=sqrt(tmp);
#else
	   w[i] = 0.5-0.5*cos(2.0*PI*(i-1)/len);
#endif
   }
   return w;
}

// hypergeom gain 
float HypergeomGain(float xx)
{
   int ind;
   float integer, frac;
   float x;
   static const float table[21] = {
      0.82157f, 1.02017f, 1.20461f, 1.37534f, 1.53363f, 1.68092f, 1.81865f,
      1.94811f, 2.07038f, 2.18638f, 2.29688f, 2.40255f, 2.50391f, 2.60144f,
      2.69551f, 2.78647f, 2.87458f, 2.96015f, 3.04333f, 3.12431f, 3.20326f};
	  x = 1*xx;
      integer = floor(2*x);
      ind = (int)integer;
      if (ind<0)
		 return 1.0;
      if (ind>19)
		 return 1.0*(1+0.1296/x);

      frac = 2*x-integer;
	  return 1.0*((1-frac)*table[ind] + frac*table[ind+1])/sqrt(x+.0001f); 
}
