/**************************************************************************
   elmo.c
 **************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"

#define DEBUG FALSE

/* -------------------- used only by elmo and dell -----------------------*/
double timelimit;
extern double timelimit;

/*-------------------------------------------------------------------------*/
void elmo(
          double t,
          int YIN, 
          int DOUT, 
          int index)
/*------------------------------------------------------------------------*/
{
    static int step;
    int i;
    double x, phix, pa, ra, sCILIA;


    if ((!activeEVANSstiffness) 
     && (!capacitiveCircuit)
     && (!extracellularGradient)
     )
        return;

    ohcscale    = 1e09;
    timelimit   = 15;

    for (i = 0;i < length; i++) 
    {
    	x   = (double) i*cl/length;

        phix = deflect(t,i,YIN);
        pa  = apical_channel_probability(phix);
        ra  = apical_resistance(pa,i);

        /* sets cilia coupling */
        if (activeEVANSstiffness)       
            sCILIA = bundle_stiffness(YIN,pa,i,x);

        if (capacitiveCircuit)            
             capacitivemembrane(ra,x,YIN,DOUT,t,i);
        else if (extracellularGradient)            
             resistivemembrane(ra,x,YIN,DOUT,t,i);

        /* show it */
        if ((!ranging)
         && (YIN == 0) 
         && (i == peakindex)
         && (t <= timelimit)) 
        {
            if (step%tds.iskip == 0) 
            {
                if ((activeEVANSstiffness)) 
                {
                    writeOHCfile(t, sCILIA,"ciliaS",3);
                }

                if ((capacitiveCircuit)||(extracellularGradient)) 
                {
                    writeOHCfile(t, Y[YIN][0][ELMO].loc[i],"V",4);
                    writeOHCfile(t, D[DOUT][0][ELMO].loc[i],"dV",5);
                }
                writeOHCfile(t, phix*1e-09,"phix",0);
                writeOHCfile(t, pa,"pa",1);
                writeOHCfile(t, ra,"ra",2);
            }
            step++;
        }
        else if (t>timelimit) 
        {
	        closeOHCfiles();
	        if (DEBUG) exit(0);
        }
    }

}
/*------------------------------------------------------------------------*/
double deflect(
               double t,
               int i,
               int YIN)
/*------------------------------------------------------------------------*/
{
    double x;
    double upper, lower;

    if (object[ACRL])
       lower = Y[YIN][0][ACRL].loc[i]; 
    else 
       lower = Y[YIN][0][BAM].loc[i];  

    if (object[TM]) 
	    upper = Y[YIN][0][TM].loc[i];
    else
	    upper = 0.0;
    
    if (MULTI)
        upper = fabs(Y[YIN][0][TM].loc[i]); 

    /* scale: convert m to nm, use height change, apply lever gain */
    /* inverted (upper - lower); correct (lower - upper) */
    x = (lower - upper) * (ohcscale * gain[i] * (stereo_h0/stereo_h[i]));

    return(x);
}

/*------------------------------------------------------------------------*/
double apical_channel_probability(double x)
/*------------------------------------------------------------------------*/
{
    double y,p;

    y=(-zP*(x-X0))/(OHC_k*OHC_T);

    /* p is 15% open at x = 0 */
    p=1/(1+exp(y));

    return(p);
}
/*------------------------------------------------------------------------*/
double bundle_stiffness(
                        int YIN,
                        double pa,
                        int i,
                        double x)
/*------------------------------------------------------------------------*/
{
    double s, volticS;


    if (coupled[TM][ACRL]) 
    {
        volticS =  staticS[i]*(Y[YIN][0][ELMO].loc[i]/Vb0[i]);
	    cS[TM][ACRL][i] = cS[ACRL][TM][i] = volticS;
    }
    else 
    {
        cS[TM][ACRL][i] = cS[ACRL][TM][i] = staticS[i];
    }

    return(s);
}
/*------------------------------------------------------------------------*/
double apical_resistance(
                         double pa,
                         int i)
/*------------------------------------------------------------------------*/
{
    double r, g;

    /* channel conductance */
    g = stereo_N[i] * gc * pa; 

    /* leak conductance is a constant defined in buildOHC; gleak = .2 */
    r = (1 / (g + gleak)); 

    r *= 1e03;  /* puts r in MOhms */

    return(r);
}
/*--------------------------------------------------------------------------*/
