#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"

/*-------------------------------------------------------------------------*/
void buildGamma(int chan)
/*-------------------------------------------------------------------------*/
{
    int i;
    FILE *fp;
    struct corey_complex cg;
    char filename[nchar];

    /* the first tone (freq[0]) determines viscosity */
    if ((!viscosity) || (chan == ST) || (chan == SM))
    { 
        for (i=0;i<length;i++) 
        {
	        Gamma0[chan][i]=1.0;
	        Gamma1[chan][i]=0.0;
	    }

        return;
    }

    strcpy(filename,viscousdir);
    strcat(filename,"/tabless");
    if ((fp = fopen(filename,"r")) == NULL)
	    error("in viscous.c \n","Open failed for tabless.\n");

    cg = readviscoustable(fp,freq[0]);

    fclose(fp);

    for (i=0;i<length;i++) 
    { /* doesn't change over length -- */
	    Gamma0[chan][i]=cg.re;
	    Gamma1[chan][i]=cg.im;
    }
}
/*-------------------------------------------------------------------------*/
struct corey_complex readviscoustable(
                                      FILE *fp,
                                      double infreq)
/*-------------------------------------------------------------------------*/
{
int rval, count;
float fr;
struct corey_complex x;

rval=0;
count=0;
while (rval!=(-1)) {
    rval = fscanf(fp,"%e %le %le %le %le",&fr,&x.re,&x.im,&x.ma,&x.ph);
    if (fr==infreq*1000) return(x);
    count++;
    } 
error("input paradigm not found in gamma table","see viscous.c");
return(x);
}
/*-------------------------------------------------------------------------*/
void initvgamma(int obj)
/*-------------------------------------------------------------------------*/
{
    int x;

    for (x=0;x<length;x++) 
        vgamma[obj][x]= 0;
}
/*-------------------------------------------------------------------------*/
void buildvgamma(int obj)
/*-------------------------------------------------------------------------*/
{
    int x;
    double wn[length];

    if (!viscosity)
    {
       for (x=0;x<length;x++) 
    	    vgamma[obj][x]= 0;
       return;
    }

    calcwnumber(obj,wn);

    if (!MULTI) 
       for (x=0;x<length;x++) 
    	    vgamma[obj][x]= (Gamma1[SM][x]/Gamma0[SM][x])*wn[x]; 
    else 
        for (x=0;x<length;x++) 
	        vgamma[obj][x]= (Gamma1[UC[obj]][x]/Gamma0[UC[obj]][x]
			    +Gamma1[LC[obj]][x]/Gamma0[LC[obj]][x])*wn[x]; 

}

/*-------------------------------------------------------------------------*/
void calcwnumber(
                 int obj,
                 double wn[length])
/*-------------------------------------------------------------------------*/
{
    int i,j,x, peak[length],n;
    int DEBUG = FALSE;
    double last, current, next, index_difference, unit_index;
    double wavelength[length],cover;

    /* calculate unit of length */
    unit_index=cl/length;

    /* calculate wavelengths -- i.e., number of millimeters */
    /* in a wave measured peak-to-peak */
    peak[0]=0;
    n=1;
    cover=0;
    for (x=1;x<length-1;x++) 
    {
	    last=Y[CURRENT][DISPLACEMENT][obj].loc[x-1];
	    current=Y[CURRENT][DISPLACEMENT][obj].loc[x];
	    next=Y[CURRENT][DISPLACEMENT][obj].loc[x+1];
    	if ((current>=last)&&(current>next)) 
        { 
		    peak[n]=x;
		    index_difference=(double) (peak[n]-peak[n-1]);
		    wavelength[n] =index_difference*unit_index;
		    cover+=wavelength[n];
		    n++;
		}
    }

    /* end condition */
    peak[n]=length;
    if (n==1) 
    {
        fprintf(stderr,"no wave yet... - kinetic viscous term is zero"); 
	    for (x=0;x<length;x++) wn[x]=0;
	    return;
	}
    else 
    {
	    index_difference=(double) (peak[n]-peak[n-1]);
	    wavelength[n] =index_difference*unit_index;
	    cover+=wavelength[n];
	}

    /* wavenumber -- inverse of wavelength */
    for (i=1;i<=n;i++) 
        for (j=peak[i-1];j<peak[i];j++) 
            wn[j]=1.0/wavelength[i];

    if (DEBUG) 
    {
        for (x=1;x<n;x++)
        {
            fprintf(stderr,"wavelength[%d] (%d to %d) = %e, total %e \n", 
	        x,peak[x-1],peak[x],wavelength[x],cover); 
	    }
    }
}
/*-------------------------------------------------------------------------*/
int MChecktime(time)
double time;
/*-------------------------------------------------------------------------*/
{
static double restarttime;
int i, rval;

rval=FALSE;
if (time==0.0) restarttime=0.0;
for (i=0;i<intone;i++) {
	if (time-restarttime>=1.0/freq[i]) {
		rval=TRUE;
		restarttime=time;
		}
	}
return(rval);
}
/*-------------------------------------------------------------------------*/
void checkgamma(
                int YF,
                int iO,
                double time)
/*-------------------------------------------------------------------------*/
{
    int minsignal;
    double newA[length][3];

    /* at the beginning of the signal period, at the first step of the rk4 */
    /* make sure that the A matrix does not need to be rebuilt because of  */
    /* changes in wavenumber */

    if (!viscosity) return;

    if ((YF==0)&&(fluidbound[iO])) 
    {
        minsignal=MChecktime(time);
	    if (minsignal) 
        {
	       buildvgamma(iO);
	       constructA(newA,iO);
	       buildLU(newA,iO);
	    }
    }
               
}
/*-------------------------------------------------------------------------*/
