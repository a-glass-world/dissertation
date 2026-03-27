#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"

#define errnum 0
/*-----------------------------------------------------------------------*/
void comment(int n)
/*-----------------------------------------------------------------------*/
{
    int i;
    int fundindex, maxindex;
    double maxloc;

    if (!tds.runfft) return;

    for (i=0;i<nobject;i++) 
    {
        if (objfft[i]<actmem)  
        { 
	        getfourierindex(&fundindex,i,1.0);
	        getplaceindex(fundindex,&maxindex,i);
	        maxloc=(double) maxindex*cl/fftlength;
	        fprintf(explain.fp,"\n\t The fundamental fourier index is %d \n",fundindex);
	        fprintf(explain.fp,"\n\t Peak response for object %d is at index %d,location %e \n",i,maxindex,maxloc); 
   	    }
    }
}
/*-----------------------------------------------------------------------*/
void buildresponse(
                   int n,
                   int p,
                   int c,
                   int o,
                   int np,
                   int stability,
                   int obsindex)
/*-----------------------------------------------------------------------*/
{
    static double calfreq;
    int fundindex,maxindex;
    int m,i,j;

    m = nmodel;

    /* set the indexes */
    if (stability != stable) 
    {
	    fundindex   = 0;
	    maxindex    = 0;
        /* rebuild the response top for all the inputs and set stability ... */
        for (i = 0; i< nobject; i++) 
        {
            if (objfft[i] < actmem)
            {
                for (j = 0; j< nscanmodel; j++)
                {
                    buildresponsetop(n, m, i, fundindex, obsindex, maxindex, p, o, c, np, stability);
                }
            }
        }
        return;
	}
    else if (tds.runfft) 
    {
	    getfourierindex(&fundindex, BAM, 1.0);
	    getplaceindex(fundindex, &maxindex, BAM);
	}
    else 
    {
	    /* this is click data and the observation place is set to 3.0 mm */
	    fundindex = 0;
	    maxindex = 20;
	}

    /* create the basic pieces of the response struct */
    for (i = 0; i < nobject; i++) 
    {
        if (objfft[i] < actmem) 
        { 
             buildresponsetop(n, m, i, fundindex, obsindex, maxindex, p, o, c, np, stability);
	    }
    }

    if (stability == stable) /* build the rest of the response struct */
    {
        if (calibrating) 
        {
            if (estimating) data[n].obsindex = maxindex;

            /* use the calibrating frequency later, for click fft analysis */
            calfreq = freq[0];
        }
        else 
        {
            for (i=0;i<nobject;i++) 
            {
	            if (objfft[i]<actmem) 
                { 
	                switch(data[n].dattype) 
                    {
    	                case liveclick: case furoclick: case postclick: case loudclick:
                        {
	                        buildclickfeat(n,m,i,fundindex,obsindex); 
			                break; 
                        } 
   	                    case isointensity:case ttones: case loudiso: case furoinout: 
                        { 
			                buildisofeat(n,m,i,fundindex,obsindex); 
			                break; 
                        }
                        case livetune: case posttune: case quicklivetune:
                        {
	       		            buildtunefeat(n,m,i,fundindex,obsindex); 
			                break;
                        } 
    	                case livegain: 
                        {
	       		            buildgainfeat(n,m,i,fundindex,obsindex); 
			                break; 
                        } 
	    	            case ttoned: 
                        {
			                error("ttone distortion case undefined",
			                "in buildresponse.c");
			                break; 
                        }
    	                case furoclickfft: 
                        {
			                buildclickfftfeat(n,m,i,calfreq,obsindex);
			                break; 
                        }
    	                default:
                        { 
			                error("datatype is undefined", "in defineinput()"); 
		                } 
	                } 
                } 
            } 
        }
    }
}
/*-------------------------------------------------------------------------*/
void buildclickfftfeat(
                       int n,
                       int m,
                       int i,
                       double calfreq,
                       int obsindex)
/*-------------------------------------------------------------------------*/
{
    int j, index;

    for (j = 0;j < response[n][m][i].nfeature;j += 2) 
    {
	    /* the fundamental index is component 1.0, the new component
	       is (featuredfreq/calibratedfrequence) */
	    response[n][m][i].feature[j] = data[n].feature[j];
	    getfourierindex(&index,i,(double) data[n].feature[j]/calfreq);
	    response[n][m][i].feature[j+1] = fftFreq[i][obsindex].mag[index]*1000000.0;
	}

}
/*-------------------------------------------------------------------------*/
void buildclickfeat(
                    int n,
                    int m,
                    int i,
                    int fundindex,
                    int obsindex)
/*-------------------------------------------------------------------------*/
{
double onset,decay,spindle,peakvel,numcycles;

buildtimefeat(i,&onset,&decay,&spindle,&peakvel,&numcycles,obsindex);

response[n][m][i].feature[0] = onset;
response[n][m][i].feature[1] = decay;
response[n][m][i].feature[2] = spindle;

/* comparing to data in micrometers per second to responses in mm/ms */
response[n][m][i].feature[3] = peakvel*1000000.0;
response[n][m][i].feature[4] = numcycles;

}
/*-------------------------------------------------------------------------*/
void buildtimefeat(
                   int i,
                   double *onset,
                   double *decay,
                   double *duration,
                   double *peakvel,
                   double *numcycles,
                   int obsindex)
/*-------------------------------------------------------------------------*/
{
    int j,x,index;
    double dbound;
    double period;
    double maxreal, minreal;
    double mintime, maxtime, peaktime;
    double spindle_start, spindle_end;
    double lastgasp,uldiv,lldiv,t,y;
    int started1;
    double halfcycle;

    /* you only get the fft time */
    index=objfft[i];

    /* find the peak to peak velocity */
    maxreal=0; minreal=0; dbound=0;

    /* use the location of CF from obsindex - during the calibration */
    x=obsindex;

    /* initialize */
    spindle_start=0;
    lastgasp=0;
    spindle_end=0;
    started1=FALSE;
    /* determines beginning of spindle */
    lldiv=2.0/5.0;
    /* determines beginning of response */
    uldiv=1.0/5.0;
    y=0;
    t=0;
    for (j=0;j<TIMESTEPS;j++) { 
	    y=ffIn[index][x].real[j];
    	    t=ffIn[index][x].time[j];
	    /* better be velocity if we're in estimation space!! */
	    if (y > maxreal) {
		    maxreal=y;
		    maxtime=t; }
	    if (y < minreal) {
		    minreal=y;
		    mintime=t; }
	    }

    *peakvel=maxreal-minreal;
    if (maxreal==0) 
    {
	    if (index == BAM) fprintf(stderr,"For BAM:\n");
	    if (index == TM) fprintf(stderr,"For TM:\n");
	    if (index == OHCSTEREO) fprintf(stderr,"For OHCSTEREO:\n");
	    if (index == ACRL) fprintf(stderr,"For ACRL:\n");
	    if (index == OHC) fprintf(stderr,"For OHC:\n");
	    error("Velocity is zero everywhere!","");
	}

    /* find the first & last time that the signal is within lldiv of the maximum */
    if (maxreal*maxreal>minreal*minreal) {
	    dbound=(maxreal*maxreal);
	    peaktime=maxtime; }
    else {
	    dbound=(minreal*minreal);
	    peaktime=mintime; }

    y=0;
    t=0;
    for (j=0;j<TIMESTEPS;j++) 
    { 
        y=ffIn[index][x].real[j];
        t=ffIn[index][x].time[j];
        if (y*y>=lldiv*dbound) 
        {
	        if (!started1) 
            {
	            spindle_start=t;
	            started1=TRUE;
	        }
	        spindle_end=t;
        }
    }
    *onset=spindle_start;
    *duration=spindle_end-spindle_start;

    /* find when it is almost at zero */
    y=0;
    t=0;
    for (j=0;j<TIMESTEPS;j++) 
    { 
	    t=ffIn[index][x].time[j];
	    y=ffIn[index][x].real[j];
	    /* yep, i'm including ringing and, no, there's no noise floor */
	    if (y*y>=dbound*.03*.03) { lastgasp=t; }
	}
    *decay=lastgasp-peaktime;

    /* find the number of cycles to spindle onset */
    y=0;
    t=0;
    halfcycle=0;
    for (j=0;j<TIMESTEPS;j++) 
    { 
	    t=ffIn[index][x].time[j];
	    y=ffIn[index][x].real[j];
	    if ((y*y>=lldiv*dbound)&&(t<=peaktime)) 
        {
		    /* count the cycles */
		    if ((y>ffIn[index][x].real[j+1])
		    &&  (y>ffIn[index][x].real[j-1])) { halfcycle++; }
		    if ((y<ffIn[index][x].real[j+1])
		    &&  (y<ffIn[index][x].real[j-1])) { halfcycle++; }
		}
	}
    *numcycles= halfcycle/2.0;

    /* check it */
    if (!RECIO) 
    {
	    period = 1.0/freq[0];
	    *numcycles = peaktime/period;
	}
}
/*-----------------------------------------------------------------------*/
void buildresponsetop(
                      int n,
                      int m,
                      int i,
                      int fundindex,
                      int obsindex,
                      int maxindex,
                      int ptype,
                      int obj,
                      int coef,
                      int np,
                      int stability)
/*-----------------------------------------------------------------------*/
{
    /* use the input values */
	response[n][m][i].stability   = stability;
	response[n][m][i].cal         = calibrating;
	response[n][m][i].fundindex   = fundindex;
	response[n][m][i].maxindex    = maxindex;
	response[n][m][i].obsindex    = obsindex;
	response[n][m][i].observe     = observe;
	response[n][m][i].freq        = freq[0];
	response[n][m][i].dB          = dB[0];
	response[n][m][i].inputnum    = n;
	response[n][m][i].signal      = signaltype;
	response[n][m][i].dattype     = datatype;

    if (estimating)
    {
	    response[n][m][i].est      = currentparm;
	    response[n][m][i].estcoef  = coef;
	    response[n][m][i].estparm  = ptype;
	    response[n][m][i].estobj   = obj;
	    response[n][m][i].newpm    = np;
        response[n][m][i].modelnum = m;

    }
    else
    {
	    response[n][m][i].est      = 0;
	    response[n][m][i].estcoef  = 0;
	    response[n][m][i].estparm  = 0;
	    response[n][m][i].estobj   = 0;
	    response[n][m][i].newpm    = FALSE;
        response[n][m][i].modelnum = 0;
    }

    /* use the data struct */
	response[n][m][i].nfeature  = data[n].nfeature;

	/* calculate the interval */
	if (calculable(n, m))
	{
		response[n][m][i].intval = response[n][m][i].est - response[n][m-1][i].est;
	}
	else
	{
		response[n][m][i].intval = 0.0;

	}
}
/*-----------------------------------------------------------------------*/
void buildisofeat(
                  int n,
                  int m,
                  int i,
                  int fundindex,
                  int obsindex)
/*-----------------------------------------------------------------------*/
{

    /* comparing to data in micrometers per second to responses in mm/ms */
    response[n][m][i].feature[0] = fftFreq[i][obsindex].mag[fundindex]*1000000.0;
}
/*-----------------------------------------------------------------------*/
void buildtunefeat(
                   int n,
                   int m,
                   int i,
                   int fundindex,
                   int obsindex)
/*-----------------------------------------------------------------------*/
{
    double p, p0;
    double mx, mx0;

    /* set to zero at base */
    p0 = fftFreq[i][0].phase[fundindex];
    p  = fftFreq[i][obsindex].phase[fundindex];

    /* data in radians */
    response[n][m][i].feature[1] = (p - p0) * 2 * M_PI;

    mx0 = fftFreq[i][length-1].mag[fundindex];
    mx0 = 0;
    mx  = fftFreq[i][obsindex].mag[fundindex];

    /* comparing to data in micrometers per second to responses in m/s */
    response[n][m][i].feature[0]=(mx - mx0)*1e06;
}
/*-----------------------------------------------------------------------*/
void buildgainfeat(
                   int n,
                   int m,
                   int i,
                   int fundindex,
                   int obsindex)
/*-----------------------------------------------------------------------*/
{
    double invelocity, indBSPL;
    double mx, mx0;

    /* comparing to data in micrometers per second to responses in mm/ms */
    mx0 =fftFreq[i][0].mag[fundindex];
    mx =fftFreq[i][obsindex].mag[fundindex];

    /* comparing to data in micrometers per second to responses in m/s */
    invelocity=(mx - mx0)*1e06;
    indBSPL=response[n][m][i].dB;
    response[n][m][i].feature[0] = measuregain(invelocity,indBSPL);
}
/*------------------------------------------------------------------------*/
double measuregain(
                   double velocity,
                   double dBSPL)
/*------------------------------------------------------------------------*/
{
    double up,answer,gain;

    up      = dBSPL / 20.0;
    answer  = exp10(up);
    gain    = velocity /( 0.02 * answer);

    return(gain);
}
/*-----------------------------------------------------------------------*/
void getfourierindex(
                     int *index,
                     int obj,
                     double comp)
/*-----------------------------------------------------------------------*/
{
    int k, oix;
    int quit = 0;

    *index = errnum;

    /* use obj */
    oix = objfft[obj];

    for (k = 0; k < TIMESTEPS/2; k++) 
    {
	    if (  (fftFreq[oix][0].freq[k]  >= comp) 
           && (fftFreq[oix][0].freq[k-1] < comp))
		    *index = k;
	}

    if (*index == errnum) 
    {
        fprintf(stderr,"in buildresponse.c, getfourierindex(),\nComponent not found\n");
        quit = query("Quit (q), continue (1).",0,1); 
    }

}
/*-----------------------------------------------------------------------*/
void getplaceindex(
                   int k,
                   int *index,
                   int obj)
/*-----------------------------------------------------------------------*/
{
    double max;
    int j, objindex;

    max      = 0.0;
    *index   = errnum;
    objindex = objfft[obj];

    /* use obj */
    for (j = 0; j < fftlength; j++) 
    {
	    if (fftFreq[objindex][j].mag[k] > max) 
        {
		    max     = fftFreq[objindex][j].mag[k];
		    *index  = j;	
	    }
    }

    /* 
    if (*index==errnum) 
        if (obj==BAM) 
            error("in buildresponse.c, getplaceindex()\n","observation place not found\n"); 
    */ 
}
/*-------------------------------------------------------------------------*/
