/**************************************************************************
   dell.c
 **************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"

#define DEBUG FALSE


double timelimit;
extern double timelimit;

/*------------------------------------------------------------------------*/
void dell(
          double t,
          int YIN,
          int index)
/*------------------------------------------------------------------------*/
{
    int i;
    static int step;
    double dL[length];

    if (!motileOHC) return;

    for (i = 0; i < length; i++) 
    {
	    dL[i] = ohcshape(Y[YIN][0][ELMO].loc[i],i);

        /* dL is in micrometers, for the force eqs., go to meters */  
	    if (coupledOHC) 
            setcouple(YIN,dL[i]*1e-06,t,index,i);
	    else 
            setzerocouple(YIN,i);

	    if (activeOHCstiffness)  setOHCstiffness(dL[i],i);
	    if (activeDEITERcouple)  changecoupling(dL[i],t,i,index);
    }

    /* show it */
    i = peakindex;
    if ((YIN == 0) && (t <= timelimit)) 
    {
        if (step%tds.iskip == 0) 
        {
            /* dL is in micrometers, for the force eqs., go to meters */
            writeOHCfile(t, (dL[peakindex] + ohc_l[i])*1e-06,"L",8);
            if (!ranging)
            {
                writeOHCfile(t, dL[i]*1e-06,"dL",7);
                if (coupledOHC)
                {
                    writeOHCfile(t, activeC[YIN][0][OHC].loc[i],"actC",9);
                    writeOHCfile(t, activeC[YIN][1][OHC].loc[i],"actCderive",10);
                }
                if (activeOHCstiffness)
                {
                    writeOHCfile(t, S[OHC][i],"S",11);
                }
                if (activeDEITERcouple)
                {
                    writeOHCfile(t, cS[OHC][BAM][i],"cS",12);
                }
            }
        }
        step++;
    }
    else if  (t > timelimit) 
    {
	    closeOHCfiles();
	    if (DEBUG) exit(0);
    }

}
/*------------------------------------------------------------------------*/
void setzerocouple(
                   int YIN,
                   int n)
/*------------------------------------------------------------------------*/
{
    activeC[YIN][0][OHC].loc[n] = 0.0;
    activeC[YIN][1][OHC].loc[n] = 0.0; 
}
/*------------------------------------------------------------------------*/
void setcouple(
               int YIN,
               double dL,
               double t,
               int index,
               int n)
/*------------------------------------------------------------------------*/
{
    double dlc = 0;
    double dlt = 0;
    double vlc = 0;

    static int started[rksteps][length];
    static double llc[rksteps][length];
    static double ltime[rksteps][length];
	    
    if (started[index][n]) 
    { 
	    dlc = llc[index][n]   - dL;
	    dlt = ltime[index][n] - t; 
	    vlc = dlc/dlt; 
	} 
    else 
    { 
	    vlc = 0; 
	    started[index][n] = TRUE; 
    } 

    llc[index][n]   = dL;
    ltime[index][n] = t; 

    activeC[YIN][0][OHC].loc[n] = dL; 
    activeC[YIN][1][OHC].loc[n] = vlc; 
}
/*------------------------------------------------------------------------*/
void changecoupling(
                    double dl,
                    double t,
                    int i,
                    int index)
/*------------------------------------------------------------------------*/
{
    static int started, count[length];
    static double tot_dl[length], avg_dl[length];
    static double holdcS[length], holdcR[length];

    int decouple                = FALSE;
    int j;

    double allowedcontraction   = 3.0; /* micrometers */
    double epsilon              = 1.0; /* tolerance */


    if ((!activeDEITERcouple)||(!coupled[OHC][BAM])) return;

    if (!started) 
    {
        /* store these for all of posterity */
        for (j = 0; j < length; j++)
        {
	        holdcS[j]   = cS[OHC][BAM][j];
	        holdcR[j]   = cR[OHC][BAM][j];
        }
	    started = TRUE;
	}
    /* compute the running average of dl (use only one step in rk4) */
    if (index == 0) 
    {
	    count[i]++;
	    tot_dl[i] += dl;
	    avg_dl[i]  = (double) tot_dl[i]/count[i];
	}


    /* decoupling logic */
    if (avg_dl[i] > allowedcontraction) decouple = TRUE;
    if (avg_dl[i] < allowedcontraction + epsilon) decouple = FALSE;


    /* perform decoupling */
    if (decouple) 
    { 
	    cS[OHC][BAM][i] = cS[BAM][OHC][i] = 0;
	    cR[OHC][BAM][i] = cR[BAM][OHC][i] = 0;
	}
    else 
    { 
	    cS[OHC][BAM][i] = cS[BAM][OHC][i] = holdcS[i];
	    cR[OHC][BAM][i] = cR[BAM][OHC][i] = holdcR[i];
	}



}
/*--------------------------------------------------------------------------*/
double ohcshape(
                double Vb,
                int i)
/*-------------------------------------------------------------------------*/
{
    /* defined by Dallos */
    double a0   = 0.0375;
    double b0   = 1.45; /* see also net.c, where expb0 is defined */
    double K    = 0.619;
    double M    = 80.0;
    double d0   = 0.082;

    double Boltzmann;
    double phi, dl, deltaVb;

    deltaVb = Vb - Vb0[i];
    Boltzmann = (1.0 / (1.0 + exp(a0 * deltaVb + b0)) - 1.0 / (1.0 + expb0));
    phi = d0 * Boltzmann;
    dl  = ohc_l[i] * M * K * phi;

    return(dl);
}
/*-------------------------------------------------------------------------*/
void setOHCstiffness(
                     double dL,
                     int i)
/*-------------------------------------------------------------------------*/
{
    double kdl, K0, Kactive;

    /* from dallos, 1997 */
    kdl = 1510;
    K0 = 200;

    Kactive= K0 + kdl * 1.0 * M_PI*ohc_r / (ohc_l[i]-dL);

    S[OHC][i] = Kactive;

}
