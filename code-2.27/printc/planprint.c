#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"


/*--------------------------------------------------------------------------*/
void planprint(
               FILE *fp,
               int n,
               int mtype) 
/*--------------------------------------------------------------------------*/
{

    fprintf(fp,"\n/********************************************************************/");
    fprintf(fp,"\n\t\t\t Model characteristics\n");
    fprintf(fp,"\n/********************************************************************/");
    list_topology(fp);
    list_geometry(fp);
    list_impedance(fp);
    list_coupling(fp);
    list_middleear(fp);
    list_input(fp,n,nmodel);
    if ((activeCOREYstiffness)||(motileOHC)) list_ohc(fp);
}
/*--------------------------------------------------------------------------*/
void list_topology(FILE *fp) 
/*--------------------------------------------------------------------------*/
{
int i,j;

	if (MULTI) {
		fprintf(fp,"\nThis is a multichannel model.\n");
        	fprintf(fp,"\nThe following channels are included: ");
        	for (i=0;i<nchannel;i++) {
			if (channel[i]) {
				fprintf(fp,"\n\t\t");
				namechannel(i,fp);
				}
			}
		fprintf(fp,"\n");
		}

	if (activeCOREYstiffness) 
		fprintf(fp,"\nThis uses Corey's active stiffness model.\n");
	if (activeEVANSstiffness) 
		fprintf(fp,"\nThis incorporates nonlinear stereocilia stiffness (from Zhang, Evans and Dallos, 1996).\n");
	if (activeOHCstiffness) 
		fprintf(fp,"\nThis includes OHC somatic stiffness.\n");
	if (motileOHC) 
		fprintf(fp,"\nThis is a motile OHC model.\n");
    if (coupledOHC)
        fprintf(fp,"\nThe OHC motility is coupled into the force-balance equations.\n");
    if ((motileOHC)&&(!coupledOHC))
        fprintf(fp,"\nThe OHC motility is not coupled into the force-balance equations.\n");
    if (activeDEITERcouple) 
		fprintf(fp,"\nThe couple is a dynamic DC couple-uncouple model.\n");
	if (capacitiveCircuit) 
		fprintf(fp,"\nThe drive to the motility is capacitive circuit model for the OHC.\n");
	if (extracellularGradient) 
		fprintf(fp,"\nThe drive to the motility is an extracellular circuit model for the OHC.\n");

        fprintf(fp,"\nThe following objects are used:");
        for (i=0;i<nobject;i++) {
		if ((object[i]) && (i != ELMO)) {
			fprintf(fp,"\n\t\t");
			nameit(i,fp);
			}
		}
	fprintf(fp,"\n");

        fprintf(fp,"\nThe following objects are analyzed using the fft: ");
        for (i=0;i<nobject;i++) {
		if (objfft[i]<actmem) {
			fprintf(fp,"\n\t\t");
			nameit(i,fp);
			}
		}
	fprintf(fp,"\n");

	if (CUTE) {
		fprintf(fp,"\nThe following objects are coupled: "); 
        for (i=0;i<nobject;i++) {
			for (j=i;j<nobject;j++) {
				if (coupled[i][j]) {
					fprintf(fp,"\n\t\t"); nameit(i,fp);
					fprintf(fp," and "); nameit(j,fp); 
                    if (i == j)
                        fprintf(fp," (adjacent couple)");
					}
				}
			}
		}
	else fprintf(fp,"\nNo objects are coupled.\n");
	fprintf(fp,"\n");
}
/*--------------------------------------------------------------------------*/
void list_geometry(FILE *fp)
/*--------------------------------------------------------------------------*/
{
    fprintf(fp,"\n\nIn the following equations, cochlear length (x) ranges from 0 to %e mm.\n",cl);

    if (!gaindefined) 
        return;

    fprintf(fp,"\nLeverage-gain -- due to ACRL rotation-change:");
    fprintf(fp,"\n\t\t%f exp(%f x)\n",gain_a11, gain_a12);

    fprintf(fp,"\nLeverage-gain -- due to stereocilia height-change:");
    fprintf(fp,"\n\t\t%f exp(%f x)\n",stereo_h0, stereo_h1);

}
/*--------------------------------------------------------------------------*/
void list_impedance(FILE *fp)
/*--------------------------------------------------------------------------*/
{
    int i;

    for (i=0;i<nobject;i++) {
	    if ((object[i]) && (i != ELMO))
        {
		    fprintf(fp,"\nThe material parameters for ");
		    nameit(i,fp);
		    fprintf(fp,": \n");
		    fprintf(fp, "\t\tmass       = %e exp(%e x)\n", mass[i][0],mass[i][1]);
		    fprintf(fp, "\t\tresistance = %e exp(%e x)\n", resist[i][0],resist[i][1]);
		    fprintf(fp, "\t\tstiffness  = %e exp(%e x)\n", stiff[i][0],stiff[i][1]);
		}
	}
}
/*--------------------------------------------------------------------------*/
void list_middleear(FILE *fp)
/*--------------------------------------------------------------------------*/
{
	if (MIDEAR) {
		fprintf(fp,"\nThe middle ear constants:\n");
		fprintf(fp,"\t\tAM         = %e\n",Am);
		fprintf(fp,"\t\tRM         = %e\n",Rm);
		fprintf(fp,"\t\tSM         = %e\n",Sm);
		fprintf(fp,"\t\tTM         = %e\n",Tm);
		fprintf(fp,"\t\tT2M        = %e\n",T2m);
		fprintf(fp,"\t\tTm2mm_ms   = %e\n",Tm2mm_ms);
		}
	else fprintf(fp,"\nNo middle ear transfer function used.\n");
}
/*--------------------------------------------------------------------------*/
void list_coupling(FILE *fp)
/*--------------------------------------------------------------------------*/
{
    int i,j;

    if (!CUTE) return;

    fprintf(fp,"\nThe coupling constants:\n");
    for (i=0;i<nobject;i++) 
    {
	    if (object[i]) 
        {
		    for (j=i;j<nobject;j++) 
            {
			    if (coupled[i][j]) 
                {
				    fprintf(fp,"\t");
				    nameit(i,fp);
				    fprintf(fp," to ");
				    nameit(j,fp);
				    fprintf(fp,":\n");
                    if ((cresist[i][j][0] != 0) || (cresist[i][j][1] != 0))
				        fprintf(fp, "\t\tresistance  = %e exp(%e x)\n", cresist[i][j][0],cresist[i][j][1]);
                    if ((cstiff[i][j][0] != 0) || (cstiff[i][j][1] != 0))
				        fprintf(fp, "\t\tstiffness   = %e exp(%e x)\n", cstiff[i][j][0],cstiff[i][j][1]);
				    fprintf(fp, "\n");
			    }
		    }
	    }
    }
}
/*--------------------------------------------------------------------------*/
void list_input(
                FILE *fp,
                int n,
                int m)
/*--------------------------------------------------------------------------*/
{
fprintf(fp, "\nThe observation place for time-evolution graphs is %e (mm) from the base.\n", peakindex*cl/length);
fprintf(fp, "\nThe input parameters for input-type %d are:\n",n);
if ((signaltype==TONE)||(signaltype==PIP)) {
	if (signaltype==TONE)
		fprintf(fp,"\t\tThis is a single tone input.\n");
	if (signaltype==PIP)
		fprintf(fp,"\t\tThis is a tone pip input.\n");
	fprintf(fp,"\t\tLevel       = %e mg/mm-sec2. \n",level[0]);
	fprintf(fp,"\t\tdB SPL      = %e (re: %e mg/mm-sec2). \n",dB[0],scale);
	fprintf(fp,"\t\tFrequency   = %e Khz",freq[0]);
	if (!estimating) {
	    fprintf(stderr,"Estimate %d, input %d is %e dB at %e Khz. \n",
	    m,n,dB[0],freq[0]);
	    }
	}
if ((signaltype==TTONED)||(signaltype==TTONES)) {
	fprintf(fp,"\t\tThis is a two-tone input.\n");
	fprintf(fp,"\t\tLevel 1     = %e mg/mm-sec2. \n",level[0]);
	fprintf(fp,"\t\tdB SPL 1    = %e (re: %e mg/mm-sec2). \n",dB[0],scale);
	fprintf(fp,"\t\tFrequency 1 = %e Khz\n",freq[0]);
	fprintf(fp,"\t\tLevel 2     = %e mg/mm-sec2. \n",level[1]);
	fprintf(fp,"\t\tdB SPL 2    = %e (re: %e mg/mm-sec2). \n",dB[1],scale);
	fprintf(fp,"\t\tFrequency 2 = %e Khz",freq[1]);
	if (!estimating) {
	    fprintf(stderr,"Estimate %d, input %d is %e dB at %e", 
	    m,n,dB[0],freq[0]);
	    fprintf(stderr,"and %e dB at %e Khz.\n",dB[1],freq[1]);
	    }
	}
if ((signaltype==CLICK)||(signaltype==CLICKFFT)) {
	fprintf(fp,"\t\tThis is a click input.\n");
	fprintf(fp,"\t\tLevel       = %e mg/mm-sec2. \n",level[0]);
	fprintf(fp,"\t\tdB SPL      = %e (re: %e mg/mm-sec2). \n",dB[0],scale);
	if (!estimating) {
	    fprintf(stderr,"Estimate %d, input %d is a %e dB CLICK. \n",
	    m,n,dB[0]);
	    }
	}
}
/*--------------------------------------------------------------------------*/
void list_ohc(FILE *fp)
/*--------------------------------------------------------------------------*/
{
    double locx;

    fprintf(fp,"\n/********************************************************************/");
    fprintf(fp,"\nOHC characteristics");
    fprintf(fp,"\n/********************************************************************/");

    fprintf(fp,"\n\nPotentials \n");
    fprintf(fp,"\t\tV0=%e mV\n",V0);
    fprintf(fp,"\t\tVm=%e mV\n\n",Vm);
    fprintf(fp,"Electrical constants \n");
    fprintf(fp,"\t\tCa0=%e pF\n",Ca0);
    fprintf(fp,"\t\tCa1=%e pF\n",Ca1);
    fprintf(fp,"\t\tCb0=%e pF\n",Cb0);
    fprintf(fp,"\t\tCb1=%e pF\n",Cb1);
    fprintf(fp,"\t\tRb0=%e MOhm\n",Rb0);
    fprintf(fp,"\t\tRb1=%e MOhm\n\n",Rb1);
    fprintf(fp,"Cilia stiffness parameters \n");
    fprintf(fp,"\t\t %e exp(%e x) N/mm\n\n",cstiff[ACRL][TM][0], cstiff[ACRL][TM][1]);
    fprintf(fp,"Number of apical channels \n");
    fprintf(fp,"\t\t%f exp(%f x)\n\n",stereo_N0, stereo_N1);
    fprintf(fp,"Cilia channel conductances \n");
    fprintf(fp,"\t\tgc = %e\n",gc);
    fprintf(fp,"\t\tgleak = %e\n\n",gleak);
    fprintf(fp,"OHC length change constants\n");
    fprintf(fp,"\t\tohc_l0 = %e \n",ohc_l0);
    fprintf(fp,"\t\tohc_l1 = %e \n\n",ohc_l1);
    fprintf(fp,"cilia resting position X0 = %e nm\n\n",X0);
    fprintf(fp,"channel open-probability force constant = %e\n",zP);

    fprintf(fp,"\n\n/********************************************************************/");
    fprintf(fp,"\nOHC characteristics at location %e (index %d)",peakindex*cl/length,peakindex);
    fprintf(fp,"\n/********************************************************************/");
    fprintf(fp,"\nElectrical constants \n");
    locx = (double) peakindex*cl/length;
    fprintf(fp,"\t\tCa[%d]=%e\n",peakindex,Ca0 * exp(Ca1*locx));
    fprintf(fp,"\t\tCb[%d]=%e\n",peakindex,Cb0 * exp(Cb1*locx));
    fprintf(fp,"\t\tRb[%d]=%e\n",peakindex,Rb0 * exp(Rb1*locx));
    fprintf(fp,"Cilia parameters \n");
    fprintf(fp,"\t\ts[%d]=%e\n",peakindex,cstiff[ACRL][TM][0]*exp(cstiff[ACRL][TM][1]*locx));
    fprintf(fp,"\t\tstereo_N[%d] = %e\n",peakindex,stereo_N[peakindex]);
    fprintf(fp,"\t\tstereo_h[%d] = %e\n",peakindex,stereo_h[peakindex]);
    fprintf(fp,"\t\tohc_l[%d] = %e \n",peakindex,ohc_l[peakindex]);
  
}
/*-------------------------------------------------------------------------*/
void nameit(
            int i,
            FILE *fp)
/*-------------------------------------------------------------------------*/
{
switch(i) {
	case BAM: {fprintf(fp,"basilar membrane "); break; }
	case TM: {fprintf(fp,"tectorial membrane "); break; }
	case OHC: {fprintf(fp,"outer hair cell "); break; }
	case OHCSTEREO: {fprintf(fp,"OHC stereocilia "); break; }
	case ACRL: {fprintf(fp,"ACRL complex "); break; }
	case ELMO: {fprintf(fp,"ELMO "); break; }
	default: {fprintf(fp,"unknown object "); break; }
	}
}
/*-------------------------------------------------------------------------*/
void namechannel(
                 int i,
                 FILE *fp)
/*-------------------------------------------------------------------------*/
{
switch(i) {
	case SM: { fprintf(fp,"scala media "); break; }
	case SS: { fprintf(fp,"spiral sulcus "); break; }
	case ST: { fprintf(fp,"scala tympani "); break; }
	}
}
/*-------------------------------------------------------------------------*/
