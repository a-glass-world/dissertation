#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gold.h>

/*-------------------------------------------------------------------------*/
void furoscan(
              int n,
              int measuretype)
/*-------------------------------------------------------------------------*/
{
static int started, hold;
int in;

if (!started) { hold=n; started=TRUE;}
in=n-hold;

switch(measuretype) {
	case furoclick: {furoclickinput(in); break;}
	case furoinout: {furoioinput(in); break;}
	case furoclickfft: {furocfftinput(in); break;}
	default: { error("measure type is not defined","furoscan()");}
	}

/* 
for (i=0;i<data[in].nfeature;i++) data[in].weight[i]=1.0;
data[in].nweight=data[in].nfeature;
weights are all equal now... */

}
/*-------------------------------------------------------------------------*/
void furoclickinput(int n)
/*-------------------------------------------------------------------------*/
{
error("input furosemide click data undefined","furoclickfft()");  
}
/*-------------------------------------------------------------------------*/
void furocfftinput(int n)
/*-------------------------------------------------------------------------*/
{

switch(n) {
/* from Ruggero and Rich, "Furosemide alters ...", J. Neuro, 1991 */
	case 0: { calibrating=TRUE; calibrate(9.8,n); break;}
	case 1: { /*  fig 5, L14, FFT of click response */
        	data[n].signal= signaltype=CLICKFFT;
        	data[n].dB= dB[0]=95;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=0.0;
		intone=0;
		observe=9.8;
        	data[n].nfeature=9;
		/* feature even gives fft component (KHz) */
		/* feature odd is velocity of fft component */
        	data[n].feature[0]=3.0; data[n].feature[1]=8.0;
        	data[n].feature[2]=5.0; data[n].feature[3]=40.0;
        	data[n].feature[4]=8.0; data[n].feature[5]=40.0;
        	data[n].feature[6]=9.0; data[n].feature[7]=30.0;
        	data[n].feature[8]=11.0; data[n].feature[9]=8.0;
		break;
		}	
	case 2: { 
        	data[n].signal= signaltype=CLICKFFT;
        	data[n].dB= dB[0]=75;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=0.0;
		intone=0;
		observe=9.8;
        	data[n].nfeature=9;
		/* feature even gives fft component (KHz) */
		/* feature odd is velocity of fft component */
        	data[n].feature[0]=3.0; data[n].feature[1]=1.0;
        	data[n].feature[2]=5.0; data[n].feature[4]=5.0;
        	data[n].feature[4]=8.0; data[n].feature[5]=5.0;
        	data[n].feature[6]=9.0; data[n].feature[7]=3.0;
        	data[n].feature[8]=11.0; data[n].feature[9]=0.5;
		break;
		} 
	default: { error("input case undefined","furocfftinput()"); }
	}
}
/*-------------------------------------------------------------------------*/
void furoioinput(int n)
/*-------------------------------------------------------------------------*/
{

switch(n) {
	/* from R&R, 1991 (ibid), fig 2, tone i/o curves, 11-28 minutes after */
	case 0: { calibrating=TRUE; calibrate(9.0,n); break;}
	case 1: {
		data[n].signal=signaltype=TONE;
		intone=1;
        	data[n].dB=dB[0]=40;
        	data[n].freq=freq[0]=9.0;
		observe=9.8;
		level[0]=scale*exp10(dB[0]/20.0);  
		data[n].nfeature=1;
		/* feature is bam velocity in micrometers per sec */
        	data[n].feature[0]=20.0; 
		break;
		}
	case 2: {
		data[n].signal=signaltype=TONE;
		intone=1;
        	data[n].dB=dB[0]=50;
        	data[n].freq=freq[0]=9.0;
		observe=9.8;
		level[0]=scale*exp10(dB[0]/20.0);  
		data[n].nfeature=1;
		/* feature is bam velocity in micrometers per sec */
        	data[n].feature[0]=80.0; 
		break;
		}
	case 3: { 
		data[n].signal=signaltype=TONE;
		intone=1;
        	data[n].dB=dB[0]=60;
        	data[n].freq=freq[0]=9.0;
		observe=9.8;
		level[0]=scale*exp10(dB[0]/20.0);  
		data[n].nfeature=1;
		/* feature is bam velocity in micrometers per sec */
        	data[n].feature[0]=250.0; 
		break;
		}
	case 4: {
		data[n].signal=signaltype=TONE;
		intone=1;
        	data[n].dB=dB[0]=70;
        	data[n].freq=freq[0]=9.0;
		observe=9.8;
		level[0]=scale*exp10(dB[0]/20.0);  
		data[n].nfeature=1;
		/* feature is bam velocity in micrometers per sec */
        	data[n].feature[0]=800.0; 
		break;
		}
	case 5: { 
        	data[n].signal=signaltype=TONE;
		intone=1;
        	data[n].dB=dB[0]=80;
        	data[n].freq=freq[0]=9.0;
		observe=9.8;
		level[0]=scale*exp10(dB[0]/20.0);  
		data[n].nfeature=1;
		/* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=1050;
        	break;
		}
	case 6: { 
        	data[n].signal=signaltype=TONE;
		intone=1;
        	data[n].dB=dB[0]=70;
        	data[n].freq=freq[0]=1.0;
		observe=9.8;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].nfeature=1;
		/* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=40;
        	break;
        	}
	case 7: { 
        	data[n].signal=signaltype=TONE;
		intone=1;
        	data[n].dB=dB[0]=80;
        	data[n].freq=freq[0]=1.0;
		observe=9.8;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].nfeature=1;
		/* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=150;
		break;
		}
	case 8: { 
        	data[n].signal= signaltype=TONE;
		intone=1;
        	data[n].freq=freq[0]=1.0;
        	data[n].dB= dB[0]=90;
		level[0]=scale*exp10(dB[0]/20.0);  
		observe=9.8;
        	data[n].nfeature=1;
		/* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=370;
		break;
		}	
	default: { error("input case undefined","furoioinput()"); }
	}
}
/*-------------------------------------------------------------------------*/

