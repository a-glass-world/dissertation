#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"		

#define incstep	0
#define decstep	1

#define step1	0
#define step2	1
#define step3	2
#define step4	3

/*----------------------------------------------------------------------*/
/* Given values for n variables Y[IN][0..n-1] find their derivatives 
   D[DERIVE1][0..n-1] at time and use the fourth-order runge-kutta method 
   to advance the solution over an interval step and return the incre-
   mented variables as Y[OUT][0..n-1], which need not be a distinct 
   array from Y[IN].  The user supplies the routines that return
   derivatives D[D_n] for Y[IN] at time -- in this case we are
   using Motility and Fluids. */
/*----------------------------------------------------------------------*/
int rk4(
        int j,
        double time,
        double step)
/* j is the object */
/*----------------------------------------------------------------------*/
{
    double instep, T2, Tstep, half, sixth;
    int x,i,stability;
    int numderive;

    half=step*0.5;
    sixth=step/6.0;
    T2=time+half;
    Tstep=time+step;

    /*initial full step calculation goes into Y[1], final calculation into Y[0] */

    if (j != ELMO) 
        numderive = nd;
    else 
        numderive = 1;

    /* derive does not change the values for Y, returns values for D */
    /* derive(j,time,instep,0,0) takes Y[0], returns D[0] */
    instep=step;
    stability=derive(j,time,instep,0,0,step1);
    if (stability==unstable) return(stability);
    for (i=0;i<numderive;i++) {
	    Y[1][i][j].stapes= Y[0][i][j].stapes + half*D[0][i][j].stapes;
	    Y[1][i][j].window= Y[0][i][j].window + half*D[0][i][j].window;
	    for (x=0;x<length;x++) 	
		    Y[1][i][j].loc[x]= Y[0][i][j].loc[x] + half*D[0][i][j].loc[x];
    	    } 

    instep=half;
    stability=derive(j,T2,instep,1,2,step2);
    if (stability==unstable) return(stability);
    for (i=0;i<numderive;i++) {
	    Y[1][i][j].stapes= Y[0][i][j].stapes + half*D[2][i][j].stapes;
	    Y[1][i][j].window= Y[0][i][j].window + half*D[2][i][j].window;
	    for (x=0;x<length;x++) 
		    Y[1][i][j].loc[x]=Y[0][i][j].loc[x] + half*D[2][i][j].loc[x];
        } 

    instep=half;
    stability=derive(j,T2,instep,1,3,step3);
    if (stability==unstable) return(stability);
    for (i=0;i<numderive;i++) {
        Y[1][i][j].stapes=Y[0][i][j].stapes + step*D[3][i][j].stapes;
        Y[1][i][j].window=Y[0][i][j].window + step*D[3][i][j].window;
        D[3][i][j].stapes += D[2][i][j].stapes;
        D[3][i][j].window += D[2][i][j].window;
        for (x=0;x<length;x++) {
	    Y[1][i][j].loc[x]= Y[0][i][j].loc[x] + step*D[3][i][j].loc[x];
	    D[3][i][j].loc[x] += D[2][i][j].loc[x];
	    } 
        } 

    stability=derive(j,Tstep,instep,1,4,step4);
    if (stability==unstable) return(stability);
    for (i=0;i<numderive;i++) {
    	    Y[0][i][j].stapes = Y[0][i][j].stapes 
	    + sixth*(D[0][i][j].stapes+D[4][i][j].stapes+2.0*D[3][i][j].stapes);
    	    Y[0][i][j].window = Y[0][i][j].window 
	    + sixth*(D[0][i][j].window+D[4][i][j].window+2.0*D[3][i][j].window);
	    for (x=0;x<length;x++) 
		    Y[0][i][j].loc[x]= Y[0][i][j].loc[x] 
		    + sixth*(D[0][i][j].loc[x]+D[4][i][j].loc[x]+2.0*D[3][i][j].loc[x]);
    	    } 

    varystep(Y[1][0][j],Y[0][0][j],time);
    varystep(Y[1][1][j],Y[0][1][j],time);
    return(stability);
} 
/*----------------------------------------------------------------------*/
int derive(
           int obj,
           double time,
           double step,
           int IN,
           int OUT,
           int indexstep)
/*----------------------------------------------------------------------*/
{
    int stability = stable;

    switch(obj) 
    {
	    case ELMO: 
        {
		    /* input Y[YIN][0][ELMO] = Vb, returns D[DOUT] = dVb/dt */
		    elmo(time,IN,OUT,indexstep);
            /* input Vb, get dL and dL/dt -- these are not integrated */
		    dell(time,IN,indexstep);
		    break;
		}
   	    case OHC: case ACRL: case BAM: case TM: 
        {
		    if (MULTI) 
		        stability = Mfluids(obj,time,IN,OUT,FALSE);
		    else 
		        stability = fluids(obj,time,IN,OUT,FALSE);
		    break;
	    }
	    default: 
        {
	 	    error("in rk.c, derive()","Unrecognized DOF\n ");
		    break;
	    }
    }
    return(stability);
}
/*----------------------------------------------------------------------*/
void varystep(
              struct motion x,
              struct motion y,
              double time)
/*----------------------------------------------------------------------*/
{
double err;
int i;
double macherr,d;
static int countok, countsaw, countdown, reduced;

if (!vary) return;
if (reduced) return; 

/* err is (estimate 1 - estimate 2)/estimate 1 */
vs.maxerr=1.0;
vs.minerr=-1.0;

macherr=.000001;

/* vary step only before fft starts, make sure the fft gets the correct stepsize */
/* thereafter reset for the next model */
if (fftstarted) {
	countok=countsaw=countdown=reduced=0;	
	return;
	}

for (i=0;i<length;i++) {
	err=(x.loc[i]-y.loc[i]);
	if (y.loc[i]!=0) {
		if (y.loc[i]<=0) d=(y.loc[i]-macherr);
		if (y.loc[i]>0) d=(y.loc[i]+macherr);
		err/=d;
		if ((err<vs.minerr)||(err>vs.maxerr)) {
		    fprintf(stderr,"Error is %e \n",err);
		    fprintf(stderr,"Changing step size from %e ",tds.timestep);
		    varyTD(decstep);
		    vs.startover=TRUE;
		    countok=0;
		    fprintf(stderr,"to %e\n",tds.timestep);
		    if (tds.timestep<tds.minstep) 
			error("\tstepsize too small\n","\tbye\n");
		    /* makes sure it never hits the minimum */
		    countdown++;
		    if (countdown>1) reduced=TRUE;
		    }
		else countok++;
	        if ((countsaw<3)&&(countok>=50*length)&&(tds.timestep<tds.normalstep)) {
		    fprintf(stderr,"Changing step size from %e ",tds.timestep);
		    varyTD(incstep);
		    fprintf(stderr,"to %e\n",tds.timestep);
		    countsaw++;
		    }
	      }
	}
vs.startover=FALSE;
} 
/*---------------------------------------------------------------------*/
void varyTD(int varytype)
/*----------------------------------------------------------------------*/
{
	double dfstable, dfskip;
	double s;

	tds.stabletime = 3.0/freq[0];
	if (varytype==incstep) tds.DOF--;
	if (varytype==decstep) tds.DOF++;
	s=128.0;
	tds.rkstepspercycle=s*tds.DOF;
	tds.iskip=(int) (tds.rkstepspercycle/tds.fftstepspercycle);
	tds.rkcyclesperstep=1.0/(s*tds.DOF);

	tds.timestep=tds.rkcyclesperstep*tds.cycletime;
	tds.fft_timestep=tds.fftcyclesperstep*tds.cycletime;
	dfstable= tds.stabletime/tds.timestep;
	tds.istable=(int) dfstable;

	if (tds.runfft) {
		dfskip=tds.iskip;
		tds.endtime = (double)tds.stabletime+(TIMESTEPS*dfskip*tds.timestep);
		}
	else    {
		tds.endtime=(double) tds.stabletime;
		}
	tds.iMAXSTEP = (int) (tds.endtime * (1.0/tds.timestep));
	fprintf(stderr,"\nEnding in %d steps, at time %f",tds.iMAXSTEP,tds.endtime);
	fprintf(stderr," (the final %d steps go thru FFT) \n", 
			TIMESTEPS*tds.iskip);
}
/*-------------------------------------------------------------------------*/
