/**************************************************************************
   buildohc.c
 **************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"

void OHCgeometry();

/*------------------------------------------------------------------------*/
void initOHC() 
/*------------------------------------------------------------------------*/
{

	int i;
	double x[length];
	int lo = 0;
	int hi = 1;

	Vm=80; 
	V0=-75; 	/* mV */

	/* apical capacitance */
	Ca0=2.2; /* 2.2 pF at base */
	Ca1=0.01; /* 2.6 pF at apex */

	/* basolateral capacitance */
	Cb0=20; /* 20 pF at base */
	Cb1=.039; /* 40 pF at apex */

	/* basolateral resistance */
	Rb0=50; /* 50 MOhm at base */
	Rb1=.06; /* 150 MOhm at apex */

	/* stereocilia stiffness constants */
	zP=290;		/* fN */

	/* the parameter index for zP is 6 */
	range[lo][6][OHC][0] =290;
	range[hi][6][OHC][0] =2900;

	X0=22;		/* deflection where channel opening is lowest (nm) */
	sl = 6.0;       /* tip-link stiffness */

	/* pivot stiffness */
	sp0 =440; 	
	sp1 =(-.08);      
	
	/* channel conductance */
	gc =.020; 	/* nS */
	gleak=.200;

	/* channel number */
	stereo_N0 = 142.0;   /* number of stereocilia per bundle */
	stereo_N1 = -.07;    /* change in number of cilia with location */
	for (i = 0; i < length; i++) {
		x[i] = ((double) i/length)*cl; 
		stereo_N[i] = (4.0/3.0)*(stereo_N0*exp(stereo_N1*x[i]));
	}

	OHCgeometry();
}
/*-------------------------------------------------------------------------*/
void OHCgeometry() 
/*------------------------------------------------------------------------*/
{
	int i;
	double x[length];

	/* geometric terms are in micrometers */

	/* terms that vary with location */
	/* stereocilia bundle height */
	stereo_h1 = 0.09;
	stereo_h0 = 1.0;

	/* OHC length */
	ohc_l0 = 20.0;
	ohc_l1 = 00.071;

	/* OHC radius */
	ohc_r = 4.0;

	/* apical membrane surface area */
	apical_a0 = 110;
	apical_a1 = 0.009;

	/* basolateral membrane surface area */
	basolateral_a0 = 2000;
	basolateral_a1 = -0.04;

	for (i = 0; i < length; i++) {
		x[i] = ((double) i/length)*cl; 
		stereo_h[i]  = stereo_h0 * exp(stereo_h1 * x[i]);
		ohc_l[i]  = ohc_l0 * exp(ohc_l1 * x[i]);
		apical_a[i] = apical_a0 * exp(apical_a1 * x[i]);
		basolateral_a[i] = basolateral_a0 * exp(basolateral_a1 * x[i]);
	}

}
