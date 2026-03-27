#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"

/*-------------------------------------------------------------------------*/
void strategyprint(FILE *fp)
/*------------------------------------------------------------------------*/
{
double dfskip, epsilon, numcycle, component_denom, component_numer;
int Fourier;

dfskip = (double) tds.iskip;

fprintf(fp,"\n\n/********************************************************************/");
fprintf(fp,"\n\t\t Solution strategy  characteristics");
fprintf(fp,"\n/********************************************************************/");
fprintf(fp,"\n\n\t Time domain in 1-D \n");
fprintf(fp,"\t Times:                         %f ms. -  %f ms. \n",
	tds.starttime,tds.endtime);
fprintf(fp,"\t Time to stability: \t\t%f\n",tds.stabletime);
fprintf(fp,"\t Steps to stability:            %d \n", tds.istable);
fprintf(fp,"\t Maximum number of steps:\t%d\n",tds.iMAXSTEP);
fprintf(fp,"\t Runge Kutta stepsize: \t\t%f ms.\n",tds.timestep);
if (tds.DOF!=1) fprintf(fp,"\t The stepsize is smaller than usual.\n");
fprintf(fp,"\t Runge Kutta steps per cycle:   %f\n",tds.rkstepspercycle);
numcycle=(double) tds.iMAXSTEP/tds.rkstepspercycle;
fprintf(fp,"\t Runge Kutta cycles:\t\t%f\n",numcycle);
if (vary) {
	fprintf(fp,"\t Variable stepper in use. \n");
	fprintf(fp,"\t Minimum stepsize:              %f ms. \n\n", tds.minstep);
	}

if (tds.runfft) {
	fprintf(fp,"\t FFT start time: \t\t %f ms.\n",tds.stabletime);
	fprintf(fp,"\t FFT steps:\t\t\t %d\n",TIMESTEPS); 
	fprintf(fp,"\t FFT steps per cycle:\t\t %f\n",tds.fftstepspercycle);
	numcycle=(double) TIMESTEPS/tds.fftstepspercycle;
	fprintf(fp,"\t FFT cycles:\t\t\t %f\n",numcycle);
	fprintf(fp,"\t FFT length:\t\t\t %d\n",fftlength);
	fprintf(fp,"\t FFT stepsize:\t\t\t %f ms.\n",tds.fft_timestep);
	fprintf(fp,"\t FFT skipstep:\t\t\t %d \n",tds.iskip);
	fprintf(fp,"\n\t All Fourier components are printed.\n");
	component_denom=TIMESTEPS*tds.fft_timestep*freq[0];
	component_numer=(TIMESTEPS/2.0)-1.0;
        epsilon=.01;  /* kicks (integer) Fourier over the edge */
	Fourier=(int) ((component_numer/component_denom)+epsilon);
	fprintf(fp,"\t Fourier components range from 0 to %d.\n",Fourier);
	if (tds.print3Dfft) fprintf(fp,"\t The 3D fft will print.\n");
	switch(tds.outputType) {
		case DISPLACEMENT: {
			fprintf(fp,"\t FFT output gives displacements.");
			break;
                       	}
       		case VELOCITY: {
               		fprintf(fp,"\t FFT output gives velocities.");
               		break;
               		}
       		default: {}
       		} 
	fprintf(fp,"\n");
	}
else fprintf(fp,"\nThe fft will not be run\n");

}
/*------------------------------------------------------------------------*/
