/*------------------------------------------------------------------------*/
/* This is the driver and 4th order runge-kutta algorithm.   		  
   Integrate starting values

	Y[step][0][object][x]
	Y[step][1][object][x]

from tds.starttime to tds.endtime.  					  */
/*------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"

/*-------------------------------------------------------------------------*/
int solve(
          int newoutput,
          int nm)
/*------------------------------------------------------------------------*/
{
    double time = tds.starttime;
    double step = tds.timestep;
    int done    = FALSE;
    int count   = 0;
    int i, stability;

    initConditions();
    stability = stable;
    do 	
    {
	    count++;
	    for (i = 0;i < nobject;i++) 
        {
            if (object[i]) 
            {
	    	    stability = rk4(i,time,step); 
		        if (stability == unstable) 
                    return(stability);
		    }
        }

	    store(START, CURRENT);
	    setoutput(count, time, newoutput, nm);
	    newoutput = FALSE;

	    if (tds.runfft) 
            fft(count, &done, time); 

	    if (!vs.startover) time += tds.timestep;
	    else fprintf(stderr,"restarting...\n");
    }
    while (!stopTD(time, &done));

    return(stability);

}
/*------------------------------------------------------------------------*/
void setoutput(
               int count,
               double time,
               int newx,
               int nm)
/*------------------------------------------------------------------------*/
{
static int clicked, snapped, skipcount;

	if (newx) 
    {
		clicked=FALSE;
		snapped=FALSE; 
    }

	if (count%100==0)
        fprintf(stderr,".");
	if (count%10000==0) 
        fprintf(stderr,"%d",count);

	/* snap - 2d = length * displacement */
	if ((time >= snapt) && (!snapped)) 
    {
		snapprint(CURRENT,DONTCARE,DONTCARE,SNAP,nm,!snapped);
		snapped=TRUE;
	}

    /* give the evolution over time - clickv and clickx */
	if (!clicked) skipcount=0; 

	if ((time < (tds.stabletime + 3.0)) && (skipcount%tds.iskip == 0)) 
    {
        snapprint(CURRENT,time,peakindex,CLICK,nm,!clicked);
		clicked=TRUE; 
    }  
    skipcount++;

}
/*------------------------------------------------------------------------*/
int stopTD(
           double time,
           int *done)
/*------------------------------------------------------------------------*/
{
	if ((time - tds.endtime) >= 0.0) 
        *done = TRUE;  

	return(*done);
}
/*------------------------------------------------------------------------*/
void store(
           int IN,
           int OUT)
/*----------------------------------------------------------------------*/
{
    int i,j,x;

    for (i=0;i<nobject;i++) 
    { 
        if (object[i]) 
        { 
            for (j=0;j<nd;j++) 
            {
	            Y[OUT][j][i].stapes = Y[IN][j][i].stapes;
	            Y[OUT][j][i].window = Y[IN][j][i].window;
	            D[OUT][j][i].stapes = D[IN][j][i].stapes;
	            D[OUT][j][i].window = D[IN][j][i].window;

	            for (x=0;x<length;x++) 
                {
		            Y[OUT][j][i].loc[x] = Y[IN][j][i].loc[x];
		            D[OUT][j][i].loc[x] = D[IN][j][i].loc[x];
	            }
            }
        }
    }

    P[OUT].stapes = P[IN].stapes;
    P[OUT].window = P[IN].window;
    for (x=0;x<length;x++) 
    {
            P[OUT].loc[x]= P[IN].loc[x];
	}
}
/*-----------------------------------------------------------------------*/
void initConditions()
/*-----------------------------------------------------------------------*/
{

    int i,j,x;

    for (i=0;i<nobject;i++) 
    { 
	    if (object[i]) 
        {
            if (i == ELMO) 
                initELMOConditions();
            else 
                initZeroConditions(i);
        }
    }

    /* needed for fluid solutions */
    for (j=0;j<rksteps;j++) 
    {
        for (x=0;x<length;x++) 
        {
            P[j].loc[x]= 0.0;
        }
    }

}
/*-----------------------------------------------------------------------*/
void initELMOConditions()
/*-----------------------------------------------------------------------*/
{
    int j, k, x;

    double pa;

    for (j = 0; j < rksteps; j++) 
    { 
        for (k = 0; k < nd; k++) 
        {
		    Y[j][k][ELMO].stapes = Y[j][k][ELMO].window = 0.0;
		    D[j][k][ELMO].stapes = D[j][k][ELMO].window = 0.0;
		    
            for (x=0;x<length;x++) 
            {
                /* Ra at 15% open probability */
                pa  = apical_channel_probability(0.0); 
                restingRa[x] = apical_resistance(pa,x);

                /* define the shape factors */
                shape_factor_a[x] = meshRb[x]/restingRa[x];

                /* init state space */
			    Vb0[x] = Y[j][k][ELMO].loc[x] = VbRest(restingRa[x],meshRb[x]); 
	    	    D[j][k][ELMO].loc[x] = 0.0;
            }
        }
    }
}
/*-----------------------------------------------------------------------*/
void initZeroConditions(int i)
/*-----------------------------------------------------------------------*/
{
    int j, k, x;

    for (j=0;j<rksteps;j++) 
    { 
        for (k=0;k<nd;k++) 
        {
		    Y[j][k][i].stapes= Y[j][k][i].window= 0.0;
		    D[j][k][i].stapes= D[j][k][i].window= 0.0;
		    
            for (x=0;x<length;x++) 
            {
			    Y[j][k][i].loc[x]= 0.0;
	    	    D[j][k][i].loc[x]= 0.0;
		    }
        }
    }
}
/*-----------------------------------------------------------------------*/
