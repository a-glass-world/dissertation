#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gold.h>

/*------------------------------------------------------------------------*/
void coupling(int cp)
/*------------------------------------------------------------------------*/
{

int lo = 0;
int hi = 1;

if (!CUTE) return;

switch(cp) {
	case neelyCP: {
		cresist[BAM][TM][0]= cresist[TM][BAM][0] =1.0;
		cresist[BAM][TM][1]= cresist[TM][BAM][1] =(-.20);
		cstiff[BAM][TM][0]= cstiff[TM][BAM][0] =1500;
		cstiff[BAM][TM][1]= cstiff[TM][BAM][1] =(-.47);
		break;
		}
	case thesis2DOFCP:  {
		cresist[BAM][TM][0]= cresist[TM][BAM][0] =1.0;
		cresist[BAM][TM][1]= cresist[TM][BAM][1] =(-.20);
		cstiff[BAM][TM][0]= cstiff[TM][BAM][0] =1500;
		cstiff[BAM][TM][1]= cstiff[TM][BAM][1] =(-.47);
		break;
		}
	case thesis2DOFactiveCP:  {
		/* deiter resistance  + aCrl stiffness coupling */
		cresist[BAM][OHC][0]= cresist[OHC][BAM][0] =5;
		cresist[BAM][OHC][1]= cresist[OHC][BAM][1] =0;
		cstiff[BAM][OHC][0]= cstiff[OHC][BAM][0] =stiff[BAM][0]/10;
		cstiff[BAM][OHC][1]= cstiff[OHC][BAM][1] =stiff[BAM][1];
		
		break;
		}
	case thesis3DOFCP:  {
		/* cilia coupling and fluid resistance */
		cresist[BAM][TM][0]= cresist[TM][BAM][0] =1.0;
		cresist[BAM][TM][1]= cresist[TM][BAM][1] =(-.20);

		cstiff[BAM][TM][0]= cstiff[TM][BAM][0] =1280;
		cstiff[BAM][TM][1]= cstiff[TM][BAM][1] =(-.07);
		break;
		}
	case thesis4DOFCP:  {

		/* fluid coupling - not sure of this term */
		cresist[ACRL][TM][0]= cresist[TM][ACRL][0] = 1.0;
		cresist[ACRL][TM][1]= cresist[TM][ACRL][1] = 0.0;
        range[lo][cACRLTMresist][nobject][0] = 1.0;
        range[hi][cACRLTMresist][nobject][0] = 10.0;
        range[lo][cACRLTMresist][nobject][1] = -.20;
        range[hi][cACRLTMresist][nobject][1] = 0.20;

		/* cilia coupling */
		cstiff[ACRL][TM][0]= cstiff[TM][ACRL][0] = 7000.0;
		cstiff[ACRL][TM][1]= cstiff[TM][ACRL][1] = 0.0;
        range[lo][cACRLTMstiff][nobject][0] = 8000.0;
        range[hi][cACRLTMstiff][nobject][0] = 30000.0;
        range[lo][cACRLTMstiff][nobject][1] = -.70;
        range[hi][cACRLTMstiff][nobject][1] = -.50;

		/* deiter stiffness  == BAM stiffness + OHC stiffness, sorta */
        cstiff[BAM][OHC][0]= cstiff[OHC][BAM][0] = 12400;
		cstiff[BAM][OHC][1]= cstiff[OHC][BAM][1] = (-0.64);
        range[lo][cOHCBAMstiff][nobject][0] = 6000;
        range[hi][cOHCBAMstiff][nobject][0] = 18000.0;
        range[lo][cOHCBAMstiff][nobject][1] = -0.3;
        range[hi][cOHCBAMstiff][nobject][1] = 0.0;

        /* deiter resistance */
		cresist[BAM][OHC][0]= cresist[OHC][BAM][0] = 50.0;
		cresist[BAM][OHC][1]= cresist[OHC][BAM][1] = 0.0;
        range[lo][cOHCBAMresist][nobject][0] = 25;
        range[hi][cOHCBAMresist][nobject][0] = 75;
        range[lo][cOHCBAMresist][nobject][1] = -0.1;
        range[hi][cOHCBAMresist][nobject][1] =  0.1;

		/* cuticular coupling */
		cstiff[ACRL][OHC][0]= cstiff[OHC][ACRL][0] = 100000;
		cstiff[ACRL][OHC][1]= cstiff[OHC][ACRL][1] = 0;
        range[lo][cACRLOHCstiff][nobject][0] = 10000;
        range[hi][cACRLOHCstiff][nobject][0] = 150000;
        range[lo][cACRLOHCstiff][nobject][1] = -0.1;
        range[hi][cACRLOHCstiff][nobject][1] = 0.1;

		/* pillar foot coupling */
		cstiff[ACRL][BAM][0]= cstiff[BAM][ACRL][0] = 10000;
		cstiff[ACRL][BAM][1]= cstiff[BAM][ACRL][1] = -0.6;
        range[lo][cACRLBAMstiff][nobject][0] = 1000;
        range[hi][cACRLBAMstiff][nobject][0] = 20000;
        range[lo][cACRLBAMstiff][nobject][1] = -0.7;
        range[hi][cACRLBAMstiff][nobject][1] = -0.1;

		/* adjacent cell coupling - phalangeal stiffness */
		cstiff[OHC][OHC][0]= 5000;
		cstiff[OHC][OHC][1]= 00;
        range[lo][cOHCOHCstiff][nobject][0] = 500;
        range[hi][cOHCOHCstiff][nobject][0] = 10000;
        range[lo][cOHCOHCstiff][nobject][1] = -0.3;
        range[hi][cOHCOHCstiff][nobject][1] = 0.0;

		/* adjacent coupling for reticular lamina - huh? */
		cstiff[ACRL][ACRL][0]= 180;
		cstiff[ACRL][ACRL][1]= -0.125;
        range[lo][cACRLACRLstiff][nobject][0] = 20;
        range[hi][cACRLACRLstiff][nobject][0] = 2000;
        range[lo][cACRLACRLstiff][nobject][1] = -0.3;
        range[hi][cACRLACRLstiff][nobject][1] = 0.0;

	    break;
		}
	default: {
		error("Couple undefined.","  Error in build/couple.c");
		break;
		}
	}
}
/*------------------------------------------------------------------------*/
