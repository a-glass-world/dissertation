#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"
/*-------------------------------------------------------------------------*/
void initTD()
/*-------------------------------------------------------------------------*/
{

	tds.minstep          = 0.000001;
	if (signaltype != CLICK )
        tds.runfft       = TRUE;
    else
        tds.runfft       = FALSE;
	tds.outputType       = VELOCITY;			
    tds.print3Dfft       = TRUE;
	tds.fftstepspercycle = 32.0;
	tds.fftcyclesperstep = 1.0/tds.fftstepspercycle;
	resetTD();
}
/*-------------------------------------------------------------------------*/
void resetTD() 		
/*-------------------------------------------------------------------------*/
{
    int i;
    double dfstable, dfskip; 

	tds.DOF = 0.0;
	for (i = 0;i < nobject;i++) 
        if ((object[i]) && (i != ELMO)) 
            tds.DOF ++;

    if (activeOHCstiffness) tds.DOF++;
    if (neelymodel) tds.DOF += 2;
    if (capacitiveCircuit)  tds.DOF += 2;

    switch(signaltype) 
    {
	    case CLICK: case CLICKFFT: case PIP:
        { 
			tds.stabletime = 0.0; 
            snapt = 2 * clickdur;
            break; 
        }
        default:
        {
            tds.stabletime = snapt; 
            break; 
        }
    }

    switch(signaltype) 
    {
		case CLICK: 
        case CLICKFFT: 
        {
			/* same as for an input of 10 KHz */
			tds.cycletime = 0.10;
			break;
        }
		case TTONES: 
        case TTONED: 
        case TONE: 
        case PIP: 
        {
			tds.cycletime = 1.0/freq[0];
		    break;
        }
	}	
	
    /* things go awry when the cycletime is 1.0 or more */
    if (tds.cycletime >= 1.0)
        tds.DOF += (int) 2*(tds.cycletime);

	tds.rkstepspercycle = 128.0 * tds.DOF;
	tds.iskip           = (int) (tds.rkstepspercycle/tds.fftstepspercycle);
	tds.rkcyclesperstep = 1.0/tds.rkstepspercycle;


	tds.timestep    = tds.rkcyclesperstep * tds.cycletime;
	tds.normalstep  = tds.rkcyclesperstep * tds.cycletime;
	tds.fft_timestep= tds.fftcyclesperstep* tds.cycletime;
	tds.starttime   = 0.0;
	dfstable        = tds.stabletime/tds.timestep;
	tds.istable     = (int) dfstable;

	if (tds.runfft) 
    {
	    dfskip      = tds.iskip;
	    tds.endtime = (double) tds.stabletime + (TIMESTEPS * dfskip * tds.timestep);

	}
	else	
    {
	    tds.endtime = (double) tds.stabletime;
	}

	if (signaltype == CLICK)
	    tds.endtime += clickdur;

	tds.iMAXSTEP    = (int) (tds.endtime * (1.0/tds.timestep));

    /* tell it */
	fprintf(stderr,"\nThis will take %d steps to time %f",tds.iMAXSTEP, tds.endtime);
	fprintf(stderr," (the last %d steps are FFTed) \n", TIMESTEPS * tds.iskip);
	fprintf(stderr,"For %f DOF, the stepsize is %f \n", tds.DOF, tds.timestep);
}
/*-------------------------------------------------------------------------*/

