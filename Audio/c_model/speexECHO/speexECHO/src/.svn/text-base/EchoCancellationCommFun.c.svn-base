
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "EchoCancellationCommFun.h"

#define PI 3.14159265358979323846



/* FFT IFFT */
static compx compxMult(compx a, compx b)
{
	compx c;
	c.real = a.real*b.real-a.imag*b.imag;
	c.imag = a.real*b.imag+a.imag*b.real;
	return c;
}

static compx compxAdd(compx a, compx b)
{
	compx c;
	c.real = a.real+b.real;
	c.imag = a.imag+b.imag;
	return c;
}

static compx compxSub(compx a, compx b)
{
	compx c;
	c.real = a.real-b.real;
	c.imag = a.imag-b.imag;
	return c;
}

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
	int LE, LE1;
	compx u, w, t;

	compx *xTemp = (compx *)malloc(N*sizeof(compx));

	int iCnt1, iCnt2, iCnt3;
    for(iCnt1=0; iCnt1<N; iCnt1++)
	{
		xTemp[iCnt1].real = x[iCnt1];
		xTemp[iCnt1].imag = 0;
	}


	RaderRank(xTemp, N);

    //for (iCnt1=0;iCnt1<N;iCnt1++)
    //{
	   //printf("Idx : %d ;%f and %f\n",iCnt1,xTemp[iCnt1].real,xTemp[iCnt1].imag);
    //}

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

   //for (iCnt1=0;iCnt1<N;iCnt1++)
   //{
	  // printf("Idx : %d ;%f and %f and %f\n",iCnt1,xTemp[iCnt1].real,xTemp[iCnt1].imag,y[iCnt1]);
   //}


	free(xTemp);

} 

void IFFT(float *x, float *y, int N, int L)
{

	int LE, LE1;
	compx u, w, t;
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
