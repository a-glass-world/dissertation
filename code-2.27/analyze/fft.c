/*-------------------------------------------------------------------------*/
/* This file runs the fast fourier transform.  */
/*-------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"


/* -----------------------------  comments  ----------------------------- */
/* fftlength is smaller, usually, than length because length must be kept 
   high for discretizing spatial mesh but once we are at the point of 
   finding the fft we only need enough discretization for the appearance 
   of smoothness -- it doesn't effect the accuracy of the fft. We decrease to 
   fftlength 64 -- anything larger and we tend to run out of memory.  */

/* FFTguts returns several frequency components.  The frequencies are 

		f   =	     i 
	     i     ------------------       i = 0, 1, ..., (Timesteps/2)-1
		       TIMESTEPS*stepsize
	
    f is in units of KHz.  */
/*-----------------------------------------------------------------------*/
void fft(
         int iCount,
         int *bDone,
         double inTime)
/*------------------------------------------------------------------------*/
{
static int iFFTCount, skipcount;
static double roll[TIMESTEPS][minmem];
int i,j;

/* initialize */
if ((inTime>=tds.stabletime)&&(!fftstarted)) 
{ 
    for (i=0;i<TIMESTEPS;i++) 
        for (j=0;j<minmem;j++) 
            roll[i][j] = 0;

	iFFTCount   = 0;
	skipcount   = 0;
	fftstarted  = TRUE;
}

if ((inTime >= tds.stabletime) && (fftstarted)) 
{ 
	if (skipcount%tds.iskip == 0) 
    {
		ConstructFFTinput(iFFTCount, CURRENT, inTime); 
		skipcount = 0;
		iFFTCount++;
	}
	skipcount++;
}

if ((iFFTCount == TIMESTEPS) && (fftstarted)) 
{
	runFFT(roll);
	*bDone = TRUE;
}

/* reset */
if (*bDone) 
    fftstarted=FALSE;
}
/*-------------------------------------------------------------------------*/
void ConstructFFTinput(
                       int iCount,
                       int YIN,
                       double inTime)     /* construct the input vector */
/*-------------------------------------------------------------------------*/
{
    int x, i1, i2, i, step;
    int countzero;

    double time, B, T;

    if (iCount >= TIMESTEPS) 
	    error("\tin fft.c, constructFFTinput()\n","\tOoops! iCount is too large.");

    time    = (double) iCount * tds.fft_timestep;
    step    = length/fftlength;

    countzero = 0 ;
    for (i = 0; i < nobject; i++) /* notice, not checking if its an integrated object */
    {
        /* just if we want an fft of it */
        i1 = objfft[i];
        if (i1 < actmem) /* actmem is the memory limitation */
        { 
	        for (x = 0; x < length; x += step) 
            {
	            i2 = x/step;
	            if ((gaindefined) && (i == OHCSTEREO)) 
                {
		            B = Y[YIN][tds.outputType][BAM].loc[x];
		            T = Y[YIN][tds.outputType][TM].loc[x];
                    if (object[OHC])
	    	            ffIn[i1][i2].real[iCount] = gain[x] * (stereo_h0/stereo_h[i])*(B - T);
                    else
	    	            ffIn[i1][i2].real[iCount] = gain[x] * (B - T);
	    	    }
	            else 
                {
	                ffIn[i1][i2].real[iCount] = Y[YIN][tds.outputType][i].loc[x];
	            }

    	        ffIn[i1][i2].im[iCount]     = 0.0;
    	        ffIn[i1][i2].time[iCount]   = inTime;
	        }
	    }
    }
}
/*-------------------------------------------------------------------------*/
void runFFT(double roll[TIMESTEPS][minmem])
/*-------------------------------------------------------------------------*/
{
    int i, x, index;
    double dat[2*TIMESTEPS+1];

    for (i=0;i<nobject;i++) 
    { 
	    index = objfft[i]; 
	    if (index < actmem) 
        { 
		    for (x = 0;x < fftlength; x++) 
            {
			    setdata(index,x,dat);
			    FFTguts(dat,TIMESTEPS,1);
			    setffFreq(index,x,dat,roll); 
		    } 
        } 
    }
} 
/*-------------------------------------------------------------------------*/
void setdata(
             int j,
             int x,
             double dat[TIMESTEPS*2 +1])
/*-------------------------------------------------------------------------*/
{
    int i;

    for (i=0;i<TIMESTEPS;i++) 
    {
	    dat[2*i+1]=ffIn[j][x].real[i];
	    dat[2*i+2]=ffIn[j][x].im[i];
    }
}
/*-------------------------------------------------------------------------*/
void FFTguts(
             double dat[2*TIMESTEPS+1],
             int nn, 
             int isign)
/*-------------------------------------------------------------------------*/
{
int  n, mmax, m, j, istep, i;
double wr, wi, wpr, wpi, wtemp, theta, dum;
double tempr, tempi;

n=2*nn;
j=1;
for (i=1;i<=n;i+=2) {
	if (j > i) {
		tempr = dat[j];
		tempi = dat[j+1];

		dat[j]  = dat[i];
		dat[j+1]= dat[i+1];
		dat[i]  = tempr;
		dat[i+1]= tempi;
		} /* end if */
	m = n/2;
	while ((m >= 2) && (j > m)) { 
		j = j-m;
		m = m/2;
		}
	j = j+m;
	} /* end ii loop */
mmax = 2;
while (n > mmax) {
	istep = 2*mmax;
	theta = 6.28318530717959/(isign*mmax);
	dum = sin(0.5*theta); 
	wpr = (-2.0*dum*dum);
	wpi = sin(theta);
	wr=1.0;
	wi=0.0;
	for (m=1;m<=mmax;m+=2) {
		for (i=m;i<=n;i+=istep) {
			j = i+mmax;
			tempr = wr*dat[j]   - wi*dat[j+1];
			tempi = wr*dat[j+1] + wi*dat[j];

			dat[j]   = dat[i] - tempr; 
			dat[j+1] = dat[i+1] - tempi;
			dat[i]   = dat[i] + tempr; 
			dat[i+1] = dat[i+1] + tempi;
			}  /* end i loop */
		wtemp = wr;
		wr = wr*wpr - wi*wpi + wr;
		wi = wi*wpr + wtemp*wpi + wi;
		} /* end m loop */
	mmax = istep;
	} /* end while */
} /* end FFT */
/*-------------------------------------------------------------------------*/
void setffFreq(
               int j,
               int x,
               double t[TIMESTEPS*2 +1],
               double roll[TIMESTEPS][minmem])
/*-------------------------------------------------------------------------*/
{
    int i,k;
    int czero;
    double divisor, multi;
    double check;
    static int started;

    if (!started) 
    {
        for (k = 0; k < iNF; k++)
        {
            FFT_max_comp[k] = 0;
            FFT_max_loc[k]  = 0;
            FFT_max_mag[k]  = 0.0;
        }
        started = TRUE;
    }

    multi = tds.fftcyclesperstep;

    for (i = 0;i<TIMESTEPS;i++) {
	    fftFreq[j][x].real[i] =  t[2*i+1]*multi;  
	    fftFreq[j][x].im[i]   =(-t[2*i+2]*multi);
	    }

    for (i=0;i<TIMESTEPS;i++) 
    {
	    fftFreq[j][x].mag[i] = sqrt(
            fftFreq[j][x].real[i] * fftFreq[j][x].real[i] +  
            fftFreq[j][x].im[i]* fftFreq[j][x].im[i]
            ); 
        if (fftFreq[j][x].mag[i] > FFT_max_mag[j]) 
        {
            FFT_max_comp[j] = i;
            FFT_max_loc[j]  = x;
            FFT_max_mag[j]  = fftFreq[j][x].mag[i];
        }
    }

    for (i=0;i<TIMESTEPS;i++) 
    {
	    check = 0;
	    if (fftFreq[j][x].real[i]!=0) 
        {
		    check = atan(fftFreq[j][x].im[i]/fftFreq[j][x].real[i]);
	    }
	    else 
        {
		    czero++;
	    }

	    if (x > 1) 
        {
	        if (((fftFreq[j][x].real[i]>=0.0)
               &&(fftFreq[j][x-1].real[i]<0.0))
               ||
               ((fftFreq[j][x].real[i]<=0.0)
               &&(fftFreq[j][x-1].real[i]>0.0))) 
            {
			    roll[i][j]++; 
		    }
	    }

	    fftFreq[j][x].phase[i] = (check/M_PI) - roll[i][j]; 
    }

    /* calculate the fourier component */
    divisor= (double) TIMESTEPS * tds.fft_timestep * freq[0];

    if ((signaltype == TONE) || (signaltype == PIP)) 
    {
        for (i = 0;i < TIMESTEPS;i++) 
        {
	        if (i < TIMESTEPS/2) 
            {
		        fftFreq[j][x].freq[i] = (double) i/divisor; 
            }
        }
    }

}
/*-------------------------------------------------------------------------*/
