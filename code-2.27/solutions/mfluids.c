
/*-------------------------------------------------------------------------*/
/* This file runs the fluid solution for the time-domain in 1D.		   */
/*-------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"

int viscosity;
extern int viscosity;

/*-------------------------------------------------------------------------*/
int Mfluids(
            int i,
            double time,
            int in,
            int out,
            int DEBUG)
/*-------------------------------------------------------------------------*/
{
    static int count;
    int stability;

    /* fluid globals */
    YF=in;
    DF=out;
    iO=i;

    /* make sure A matrix is OK (changes dynamically for viscous models) */
    if (viscosity) checkgamma(YF,iO,time);

    CalculateG();
    CalculateGc();
    if (fluidbound[iO]!=TRUE) 
    {
	    QuickDerive(); 
	}
    else 
    {
	    GetSignal(time);
	    totsignal=MBoundaryConditions();
	    if (MCalculatePHI())
		    MCalculateK(totsignal);
	    else
		    CalculateK(totsignal);
	    CalculateP();
	    Derive();
	    if (MIDEAR) 
        {
	        if (incident[iO][SM]) DeriveStapes();
	        if (incident[iO][ST]) MDeriveRWindow(); 
	    }
    }

    count++;
    stability=CheckUnstable(count); 

    if (DEBUG) fluiddebug(time,YF,DF,iO);

    return(stability);
}
/*-------------------------------------------------------------------------*/
double MBoundaryConditions() 
/*-------------------------------------------------------------------------*/
{
    double signal1 = 0.0;
    double signal2 = 0.0;

    signal1 = calcBC(UC[iO]);
    signal2 = calcBC(LC[iO]);

    return(signal1 + signal2);
}
/*-------------------------------------------------------------------------*/
double calcBC(int chan) 
/*-------------------------------------------------------------------------*/
{
    double fsignal = 0.0;

    switch(chan) {
        case SM: case ST: {
	        if (MIDEAR) 
            {
	            Gme = T2m*Rm*Y[YF][1][iO].stapes + Sm*Y[YF][0][iO].stapes;
	            fsignal = density*(Gme - Tm*Am*finput)/Tm2mm_ms; 
	        }
	        else 
            {
                fsignal = finput; 
            }
            break;
	    } 
        case SS:  default: {
	        fsignal=0.0;
	        break;
	    }
    }

    return(fsignal);
}
/*-------------------------------------------------------------------------*/
int MCalculatePHI()
/*-------------------------------------------------------------------------*/
{
    int i, phi_calculated;

    phi_calculated = FALSE;

    for (i=0;i<length;i++) 
        PHI[i]=0.0;

    /* phi is the sum of accelerations of membranes also incident to channels
       bordering the current object. */
    for (i = 0;i < nobject;i++) 
    {
        if  ((i != iO) && (object[i])) 
        {
            if (incident[i][UC[iO]]) 
                doPHI(UC[iO],i);
            if (incident[i][LC[iO]]) 
                doPHI(LC[iO],i);

            phi_calculated = TRUE;
        }
	}

    return(phi_calculated);
}
/*-------------------------------------------------------------------------*/
void doPHI(
           int chan,
           int obj)
/*-------------------------------------------------------------------------*/
{
    int k;

    for (k=0;k<length;k++) 
    {
    /* the acceleration of membrane is last current derivative */
	    PHI[k] += (density*beta[obj][k])/(Gamma0[chan][k]*area[chan][k])
		    *D[CURRENT][1][obj].loc[k];
	}
}
/*-------------------------------------------------------------------------*/
void MCalculateK(double insignal)
/*-------------------------------------------------------------------------*/
{
    int i;

    /* the first term of K includes the boundary condition at the base */
    K[0] = 0.5*delta[0]*(alpha[iO][0]*(G[0]+Gc[0]) - PHI[0]) + insignal;
    for (i=1; i<length; i++) 
	    K[i]=0.5 *(mesh[i+1]-mesh[i-1])*(alpha[iO][i]*(G[i]+Gc[i])-PHI[i]);

}
/*-------------------------------------------------------------------------*/
void MDeriveRWindow()
/*-------------------------------------------------------------------------*/
{
    D[DF][1][iO].window= (area[ST][0]*P[DF].loc[0]-Gwindow)/Tm2mm_ms;
    D[DF][0][iO].window=Y[YF][1][iO].window;
}
/*-------------------------------------------------------------------------*/
