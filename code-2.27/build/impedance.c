#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"

int lo = 0;
int hi = 1;

/*----------------------------------------------------------------------------*/
void impedance(int Ztype)
/*----------------------------------------------------------------------------*/
{

    /* set range in KHz - for chinchilla */
    hifreq=25.0;
    lofreq=.09;
    switch(Ztype) 
    {
	    case viergeverZ: 
        {   /* units are mm mg ms */
		    /* used to compare this solution strategy
		    to Viergever's complete LG solution - see his
		    thesis, Figures 5.2.4
		    values  
                R0= 13 uNs/mm3
			    M0= 1.5 mg/mm2
			    S0= 10 N/mm3
			    R1= 0 mm-1
			    M1= 0 mm-1
			    S1=-0.2 mm-1
		    density is 1 mg/mm3
		    cochlear length is 35 mm (human)
		    channel height is 1 mm 
		    */
		    resist[BAM][0]=13.0;
		    resist[BAM][1]=.00;
		    mass[BAM][0]= 1.5;
		    mass[BAM][1]= .00;
		    stiff[BAM][0]= 10000.0;
		    stiff[BAM][1]= (-0.2);
            if (motileOHC) initOHC();
		    break;
		}
	    case diepenZ: { /* units are mm mg ms */
		    resist[BAM][0]=5.0;
		    resist[BAM][1]=(-0.15);
		    mass[BAM][0]= .5;
		    mass[BAM][1]= .00;
		    stiff[BAM][0]= 20000.0;
		    stiff[BAM][1]= (-0.3);
		    break;
		}
	    case thesisZ1DOF: { 
		    CPimpedance();
            if (motileOHC) initOHC();
		    break;
		}
	    case thesisZ2DOF:  {
		    CPimpedance();
		    TMimpedance();
		    break;
		}
	    case thesisZ2DOFactive:  {
		    CPimpedance();
		    OHCimpedance();
		    initOHC();
		    break;
		}	
	    case thesisZ3DOF:  {
		    CPimpedance();
		    TMimpedance();
		    break;
		}	
	    case thesisZ4DOF:  
        {
            /* this is the lever gain */
            gain_a11 = 1.2;
            gain_a12 = 0.06;
		    gaindefined = TRUE;

		    BAMimpedance();
		    TMimpedance();
		    ACRLimpedance();
		    OHCimpedance();
            initOHC();
		    if (object[ELMO]) ELMOimpedance();
		    break;
		}
	    case neelyZ_PSA2DOF: {
		    resist[BAM][0]= .50;
		    resist[BAM][1]= .05;
		    mass[BAM][0]= .045;
		    mass[BAM][1]= .08;
		    stiff[BAM][0]= 4600.0;
		    stiff[BAM][1]= (-.32);

 		    resist[TM][0] = 9.0e-01;
        	resist[TM][1] = (-3.0e-01);
        	stiff[TM][0]  = 3.8e+02;
        	stiff[TM][1]  = (-5.26e-01);
        	mass[TM][0]   = 1.8e-02;
        	mass[TM][1]   = 3.0e-02;

		    gain_a11=1.2;
    		gain_a12=(-.05);
		    gaindefined = TRUE;
            break;
		}
	    default: {
		    error("in impedance.c \n","undefined case \n");
		}
	}
}
/*------------------------------------------------------------------------*/
void CPimpedance()
/*------------------------------------------------------------------------*/
{

	/* start with plate values */
	mass[BAM][0]= 0.64;
	mass[BAM][1]= 0.02;
	stiff[BAM][0]= 10000.0;
	stiff[BAM][1]= (-0.64);
	resist[BAM][0]=25.0;
	resist[BAM][1]=-.275;

	range[lo][MASS][BAM][0] = .5;
	range[hi][MASS][BAM][0] = .8;
	range[lo][MASS][BAM][1] = .00;
	range[hi][MASS][BAM][1] = .00;

	range[lo][STIFF][BAM][0] = 5000.0;
	range[hi][STIFF][BAM][0] = 20000.0;
	range[lo][STIFF][BAM][1] = -.5;
	range[hi][STIFF][BAM][1] = -.7;

	range[lo][RESIST][BAM][0] = 15.0;
	range[hi][RESIST][BAM][0] = 50.0;
	range[lo][RESIST][BAM][1] = -0.3;
	range[hi][RESIST][BAM][1] = 0.0;

}
/*------------------------------------------------------------------------*/
void BAMimpedance()
/*------------------------------------------------------------------------*/
{
   
	mass[BAM][0]= 0.9;
	mass[BAM][1]= .005;
	stiff[BAM][0]= 10000.0;
	stiff[BAM][1]= (-0.64);
	resist[BAM][0]=18.0;
	resist[BAM][1]=-.28;

	range[lo][MASS][BAM][0] = .5;
    range[hi][MASS][BAM][0] = 2.0;
	range[lo][MASS][BAM][1] = 0.0;
    range[hi][MASS][BAM][1] = 0.1;

	range[lo][STIFF][BAM][0] = 6000.0;
    range[hi][STIFF][BAM][0] = 25000.0;
	range[lo][STIFF][BAM][1] = -.7;
    range[hi][STIFF][BAM][1] = -.5;

	range[lo][RESIST][BAM][0] = 5;
    range[hi][RESIST][BAM][0] = 50;
	range[lo][RESIST][BAM][1] = -.5;
    range[hi][RESIST][BAM][1] = -.1;

}
/*------------------------------------------------------------------------*/
void TMimpedance()
/*------------------------------------------------------------------------*/
{

    mass[TM][0]   = 0.44;
    mass[TM][1]   = 0.00;
    stiff[TM][0]  = 20.0;
    stiff[TM][1]  = -0.15;
 	resist[TM][0] = 0.630;
    resist[TM][1] = 0.035;

	range[lo][MASS][TM][0] = 0.1;
    range[hi][MASS][TM][0] = 2.0;
	range[lo][MASS][TM][1] = -0.1;
    range[hi][MASS][TM][1] =  0.1;

	range[lo][STIFF][TM][0] = 2;
    range[hi][STIFF][TM][0] = 200;
	range[lo][STIFF][TM][1] = -0.3;
    range[hi][STIFF][TM][1] =  0.0;

    range[lo][RESIST][TM][0] = .20;
    range[hi][RESIST][TM][0] = 10.0;
	range[lo][RESIST][TM][1] = -0.1;
    range[hi][RESIST][TM][1] =  0.1;

}
/*------------------------------------------------------------------------*/
void ACRLimpedance()
/*------------------------------------------------------------------------*/
{

    mass[ACRL][0]   = 0.56;
    mass[ACRL][1]   = 0.005;
    stiff[ACRL][0]  = 270.0;
    stiff[ACRL][1]  = -0.125;
 	resist[ACRL][0] = 2.44;
    resist[ACRL][1] = -0.065;

	range[lo][MASS][ACRL][0] =  0.1;
    range[hi][MASS][ACRL][0] =  2.0;
	range[lo][MASS][ACRL][1] = -0.1;
    range[hi][MASS][ACRL][1] =  0.1;

	range[lo][STIFF][ACRL][0] = 20;
    range[hi][STIFF][ACRL][0] = 2000;
	range[lo][STIFF][ACRL][1] = -0.3;
    range[hi][STIFF][ACRL][1] =  0.0;

	range[lo][RESIST][ACRL][0] = 5;
    range[hi][RESIST][ACRL][0] = 50;
	range[lo][RESIST][ACRL][1] = -0.3;
    range[hi][RESIST][ACRL][1] =  0.1;

}
/*-------------------------------------------------------------------------*/
void OHCimpedance()
/*------------------------------------------------------------------------*/
{

	mass[OHC][0]= .16;
	mass[OHC][1]= .13;
	stiff[OHC][0]= 12400;
	stiff[OHC][1]= -.04;
	resist[OHC][0]= 50.0;
	resist[OHC][1]= 0.0;

	range[lo][MASS][OHC][0] = 0.05;
    range[hi][MASS][OHC][0] = 1.0;
	range[lo][MASS][OHC][1] = 0;
    range[hi][MASS][OHC][1] = .2;

	range[lo][STIFF][OHC][0] = 1240;
    range[hi][STIFF][OHC][0] = 25000;
	range[lo][STIFF][OHC][1] = -0.1;
    range[hi][STIFF][OHC][1] =  0.1;

	range[lo][RESIST][OHC][0] = 5;
    range[hi][RESIST][OHC][0] = 50;
	range[lo][RESIST][OHC][1] = -0.1;
    range[hi][RESIST][OHC][1] =  0.1;

}
/*-------------------------------------------------------------------------*/
void ELMOimpedance()
/*------------------------------------------------------------------------*/
{

	/* this is a hack to get around the finite element mesh calculated on the mechanical features */
	mass[ELMO][0]= .05;
	mass[ELMO][1]= .07;
	stiff[ELMO][0]= 12400;
	stiff[ELMO][1]= -.04;
	resist[ELMO][0]= 50.0;
	resist[ELMO][1]= 0.0;

	/* never never never estimate these phony values */
	range[lo][MASS][ELMO][0] = range[hi][MASS][ELMO][0] = mass[ELMO][0];
	range[lo][MASS][ELMO][1] = range[hi][MASS][ELMO][1] = mass[ELMO][1];
	range[lo][STIFF][ELMO][0] = range[hi][STIFF][ELMO][0] = stiff[ELMO][0];
	range[lo][STIFF][ELMO][1] = range[hi][STIFF][ELMO][1] = stiff[ELMO][1];
	range[lo][RESIST][ELMO][0] = range[hi][RESIST][ELMO][0] = resist[ELMO][0];
	range[lo][RESIST][ELMO][1] = range[hi][RESIST][ELMO][1] = resist[ELMO][1];

}

