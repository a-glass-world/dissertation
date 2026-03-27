#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"

/*-------------------------------------------------------------------------*/
void postscan(
              int n,
              int measuretype)
/*-------------------------------------------------------------------------*/
{
static int started, hold;
int in;

if (!started) { hold=n; started=TRUE;}
in=n-hold;

switch(measuretype) {
	case posttune: { posttuningcurve(in); break;}
	case postclick: {postclickinput(in); break;}
	default: { error("measure type is not defined","postscan()");}
	}
}
/*-------------------------------------------------------------------------*/
void posttuningcurve(int n)
/*-------------------------------------------------------------------------*/
{

    /* In figure 8 of the Robles et al 1986 paper one can trace the
       average values for the "passive" or damaged response.  A quick
       fit of the data gives (for the "insensitive" response) the following:
       an isovelocity response of 0.1 mm/sec for
           KHz     dB SPL   Phase
	    1       80      -1.3
	    3       68      -1.5
	    5       60      -2.5
	    6       63      -3
	    7       65      -3.8  -- not included 
	    8       67      -5.0
	    9.0     72      -7.8
	    10.0    80      -8    -- not included 

    */

    switch(n) {
	    /* the CF for that data is 5.0 KHz */
	    case 0: { 
            calibrating=TRUE; calibrate(5.0,n); break;}
	    case 1: { 
            data[n].signal=signaltype=TONE;
            data[n].dB=dB[0]=74;
		    level[0]=scale*exp10(dB[0]/20.0);  
            data[n].freq=freq[0]=1.0;
		    observe=5.0;
		    intone=1;
	        data[n].nfeature=2;
            data[n].feature[0]=100.0; 
            data[n].feature[1]=-1.3;
		    break;
		    }
	    case 2: { 
            data[n].signal=signaltype=TONE;
            data[n].dB=dB[0]=68;
		    level[0]=scale*exp10(dB[0]/20.0);  
            data[n].freq=freq[0]=3;
		    observe=5.0;
		    intone=1;
		    data[n].nfeature=2;
            data[n].feature[0]=100.0; 
            data[n].feature[1]=-1.5;
		    break;
		    }
	    case 3: { 
            data[n].signal=signaltype=TONE;
            data[n].dB=dB[0]=60;
		    level[0]=scale*exp10(dB[0]/20.0);  
            data[n].freq=freq[0]=5;
		    observe=5.0;
		    intone=1;
		    data[n].nfeature=2;
            data[n].feature[0]=100.0; 
            data[n].feature[1]=-2.5;
		    break;
		    }
	    case 4: { 
            data[n].signal=signaltype=TONE;
            data[n].dB=dB[0]=63;
		    level[0]=scale*exp10(dB[0]/20.0);  
            data[n].freq=freq[0]=6.0;
		    observe=5.0;
		    intone=1;
		    data[n].nfeature=2;
            data[n].feature[0]=100.0; 
            data[n].feature[1]=-3.0;
		    break;
		    }
	    case 5: { 
            data[n].signal=signaltype=TONE;
            data[n].dB=dB[0]=67;
		    level[0]=scale*exp10(dB[0]/20.0);  
            data[n].freq=freq[0]=8.0;
		    observe=5.0;
		    intone=1;
		    data[n].nfeature=2;
            data[n].feature[0]=100.0; 
            data[n].feature[1]=-5.0;
		    break;
		    }
	    case 6: { 
            data[n].signal=signaltype=TONE;
            data[n].dB=dB[0]=72;
		    level[0]=scale*exp10(dB[0]/20.0);  
            data[n].freq=freq[0]=9.0;
		    observe=5.0;
		    intone=1;
		    data[n].nfeature=2;
            data[n].feature[0]=100.0; 
            data[n].feature[1]=-7.8;
		    break;
		    }
	    default: { 
            error("input case undefined","posttuningcurve()"); 
        }
    }
}
/*-------------------------------------------------------------------------*/
void postclickinput(int n)
/*-------------------------------------------------------------------------*/
{

    /* the data for this input is from data file l13pmwfs.txt by Nola Rich.
       the figure showing data is results/1DOF/clickgraphs/figClickPassive.m  */

    switch(n) 
    {
    /* this is calibrated to 6.0 KHz because the place of 6.0 CF
       passive is equal to the place of 9.0 CF active */
	    case 0: { calibrating=TRUE; calibrate(6.0,n); break;}
	    case 1: { 
		    /* animal L13, time-domain clicks */
            data[n].signal=signaltype=CLICK;
            data[n].dB=dB[0]=66;
		    level[0]=scale*exp10(dB[0]/20.0);  
		    observe=6.0;
		    intone=0;
            data[n].freq=0.0;
		    data[n].nfeature=5;
		    /* feature is time in ms -- onset time */ 
            data[n].feature[0]=.4; 
		    /* feature is time in ms -- decay */ 
            data[n].feature[1]=.5; 
		    /* feature is time in ms -- duration of spindle */ 
            data[n].feature[2]=.5; 
		    /* feature is peak-to-peak velocity in micrometers/sec */ 
            data[n].feature[3]=90.0; 
		    /* feature is number of cycles to near-maximum */ 
            data[n].feature[4]=1.5; 
		    break;
		    }
	    case 2: { 
		    /* animal L13, time-domain clicks */
            data[n].signal=signaltype=CLICK;
            data[n].dB=dB[0]=76;
		    level[0]=scale*exp10(dB[0]/20.0);  
		    observe=6.0;
		    intone=0;
            data[n].freq=0.0;
		    data[n].nfeature=5;
		    /* feature is time in ms -- onset time */ 
            data[n].feature[0]=.4; 
		    /* feature is time in ms -- decay */ 
            data[n].feature[1]=.5; 
		    /* feature is time in ms -- duration of spindle */ 
            data[n].feature[2]=.5; 
		    /* feature is peak-to-peak velocity in micrometers/sec */ 
            data[n].feature[3]=290.0; 
		    /* feature is number of cycles to near-maximum */ 
            data[n].feature[4]=1.5; 
		    break;
		    }
	    case 3: { 
		    /* animal L13, time-domain clicks */
            data[n].signal=signaltype=CLICK;
            data[n].dB=dB[0]=86;
		    level[0]=scale*exp10(dB[0]/20.0);  
		    observe=6.0;
		    intone=0;
            data[n].freq=0.0;
		    data[n].nfeature=5;
		    /* feature is time in ms -- onset time */ 
            data[n].feature[0]=.4; 
		    /* feature is time in ms -- decay */ 
            data[n].feature[1]=.5; 
		    /* feature is time in ms -- duration of spindle */ 
            data[n].feature[2]=.5; 
		    /* feature is peak-to-peak velocity in micrometers/sec */ 
            data[n].feature[3]=800.0; 
		    /* feature is number of cycles to near-maximum */ 
            data[n].feature[4]=1.5; 
		    break;
		    }
	    default: { 
            error("input case undefined","postclickinput()"); 
        }
    }
}
/*-------------------------------------------------------------------------*/
