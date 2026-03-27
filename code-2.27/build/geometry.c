#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"

/*-------------------------------------------------------------------------*/
void geometry(int GEOtype) 
/*-------------------------------------------------------------------------*/
{
int ih;
double scale0;

density=1.0; 		/* rho is 1.0 mg/mm3 */
switch(GEOtype) {
	case thesisGEO: {
		/* effective height/width correction */
		scale0=.5;
		/* chinchilla */
		cl=18.8;

		/* fit to LIM 1980: 230 um < width < 370 um */
		/* using b0*exp(b1*x) */
		plate[BAM] = TRUE;
		width[BAM][0]=.230; /* mm */
		width[BAM][1]=.026;

		plate[TM] = TRUE;
		width[TM][0]=.06; /* mm */
		width[TM][1]=.026;

		plate[ACRL] = TRUE;
		width[ACRL][0]=.2; /* mm */
		width[ACRL][1]=.026;

		plate[OHC] = FALSE;
		plate[ELMO] = FALSE;
		/* Bohne data */
		height[SM][0]=scale0*.45;
		height[SM][1]=0.57;
		height[SM][2]=(-0.02);

		height[ST][0]=scale0*.33;
		height[ST][1]=-0.67;
		height[ST][2]=(-0.02);

		height[SS][0]=scale0*.1;
		height[SS][1]=0.0;
		height[SS][2]=0.0;
		break;
		}
	case diepenGEO: { 
	/* used to compare this solution strategy to 
	   Diependaal et al 1987 */
		cl=35.0;
		height[SM][0]=1.0;
		height[SM][1]=0.0;
		height[SM][2]=0.0;
		height[ST][0]=height[SM][0];
		height[ST][1]=0.0;
		height[ST][2]=0.0;
		break;
		}
	case viergeverGEO: {
	/* used to compare this solution strategy to 
	   Viergever's complete LG solution - see his 
	   thesis, Figures 5.2.2 through 5.2.8.  
	   values  
		R0= (varying) 
		M0=1.5 mg/mm2 
		S0=10 N/mm 
		R1=0 mm-1 
		M1=0 mm-1 
		S1=-0.2 mm-1 
		density is 1 mg/mm3 
		cochlear length is 35 mm (human) 
		channel height is 1 mm */
		cl=35.0;
		ih=query("Real(0) or ideal(1) geometries",0,1);
		if (ih==1) {
			height[SM][0]=1.0;
			height[SM][1]=0.0;
			height[SM][2]=0.0;
			height[ST][0]=height[SM][0];
			height[ST][1]=0.0;
			height[ST][2]=0.0;
			width[BAM][0]=1.0;
			width[BAM][1]=0.0;
			fprintf(stderr,"The next input is divided by 10.\n");
		        ih=query("Input the constant effective height",0,100);
		        height[SM][0]=((double) ih)/10.0;
		        fprintf(stderr,"The scala height is now %f \n",
		        height[SM][0]);
		        }
		else {
			width[BAM][0]=.20;
			width[BAM][1]=0.02;
			fprintf(stderr,"The next input is divided by 10.\n");
		        ih=query("Input the effective height coefficient",5,20);
		        height[SM][0]=.45*((double) ih)/10.0;
			height[ST][0]=.33*((double) ih)/10.0;
		        fprintf(stderr,"Scalae heights are %f and %f \n",
		        height[SM][0], height[ST][0]);
			height[SM][1]=0.57;
			height[SM][2]=(-0.02);
			height[ST][1]=1.5;
			height[ST][2]=(-0.02);
		        }
		break;
		}
	case neelyGEO: {
		cl=22.5;
		width[BAM][0]=1.0;
		width[BAM][1]=0.0;
		width[TM][0]=1.0;
		width[TM][1]=0.0;
		height[SM][0]=1.0; 
		height[SM][1]=0.0;
		height[SM][2]=0.0;
		height[ST][0]=height[SM][0];
		height[ST][1]=0.0;
		height[ST][2]=0.0;
		height[SS][0]=1.0;
		height[SS][1]=0.0;
		height[SS][2]=0.0;
		break;
		}
	case chap3GEO: {
	/* old values - used for chapter 3 */
		cl=22.5;
		width[BAM][0]=.20;
		width[BAM][1]=0.02;
		width[TM][0]=.06;
		width[TM][1]=0.0;

		height[SM][0]=.45;
		height[SM][1]=0.57;
		height[SM][2]=(-0.02);

		height[ST][0]=.33;
		height[ST][1]=1.5;
		height[ST][2]=(-0.02);

		height[SS][0]=.07;
		height[SS][1]=0.0;
		height[SS][2]=0.0;
		break;
		}
	default:{ }
	}

}
/*-------------------------------------------------------------------------*/
