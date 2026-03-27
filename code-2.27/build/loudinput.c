#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gold.h>

/*-------------------------------------------------------------------------*/
void loudscan(
              int n,
              int measuretype)
/*-------------------------------------------------------------------------*/
{
static int hold,started;
int in;

if (!started) { hold=n; started=TRUE;}
in=n-hold;

switch(measuretype) {
	case loudiso: {loudisoinput(in); break;}
	case loudclick: {loudclickinput(in); break;}
	default: { error("measure type is not defined","loudscan()");}
 	}	

/* 
for (i=0;i<data[in].nfeature;i++) data[in].weight[i]=1.0;
data[in].nweight=data[in].nfeature;
weights are all equal now... */

}
/*-------------------------------------------------------------------------*/
void loudisoinput(int n)
/*-------------------------------------------------------------------------*/
{
switch(n) {
	case 0: { calibrating=TRUE; calibrate(10.0,n); break;} 
	case 1: { /* observation checked -- Figure 7 in Ruggero 1995 */
		/* animal L110, Isointensity curve for 30 dB at 9 KHz */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=30;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=9.0;
		observe=10.0;
		intone=1;
		data[n].nfeature=1;
		/* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=15; 
		break;
		}
	case 2: { /* observation checked -- Figure 7 in Ruggero 1995 */
		/* animal L110, Isointensity curve for 30 dB at 10 KHz */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=30;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=10.0;
		observe=10.0;
		intone=1;
		data[n].nfeature=1;
		/* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=30; 
		break;
		}
	case 3: { /* observation checked -- Figure 7 in Ruggero 1995 */
		/* animal L110, Isointensity curve */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=30;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=11.0;
		observe=10.0;
		intone=1;
		data[n].nfeature=1;
		/* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=20; 
		break;
		}
	case 4: { /* observation checked -- Figure 7 in Ruggero 1995 */
		/* animal L110, Isointensity curve */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=50;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=7.0;
		observe=10.0;
		intone=1;
		data[n].nfeature=1;
		/* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=25; 
		break;
		}
	case 5: { /* observation checked -- Figure 7 in Ruggero 1995 */
		/* animal L110, Isointensity curve */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=50;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=8.0;
		observe=10.0;
		intone=1;
		data[n].nfeature=1;
		/* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=100; 
		break;
		}
	case 6: { /* observation checked -- Figure 7 in Ruggero 1995 */
		/* animal L110, Isointensity curve */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=50;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=9.0;
		observe=10.0;
		intone=1;
		data[n].nfeature=1;
		/* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=110; 
		break;
		}
	case 7: { /* observation checked -- Figure 7 in Ruggero 1995 */
		/* animal L110, Isointensity curve */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=50;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=10.0;
		observe=10.0;
		intone=1;
		data[n].nfeature=1;
		/* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=250; 
		break;
		}
	case 8: { /* observation checked -- Figure 7 in Ruggero 1995 */
		/* animal L110, Isointensity curve */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=50;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=11.0;
		observe=10.0;
		intone=1;
		data[n].nfeature=1;
		/* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=120; 
		break;
		}
	case 9: { /* observation checked -- Figure 7 in Ruggero 1995 */
		/* animal L110, Isointensity curve */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=80;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=9.0;
		observe=10.0;
		intone=1;
		data[n].nfeature=1;
		/* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=1000; 
		break;
		}
	case 10: { /* observation checked -- Figure 7 in Ruggero 1995 */
		/* animal L110, Isointensity curve */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=80;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=10.0;
		observe=10.0;
		intone=1;
		data[n].nfeature=1;
		/* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=1000; 
		break;
		}
	case 11: { /* observation checked -- Figure 7 in Ruggero 1995 */
		/* animal L110, Isointensity curve */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=80;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=11.0;
		observe=10.0;
		intone=1;
		data[n].nfeature=1;
		/* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=400; 
		break;
		}
	default: { error("input case undefined","loudisoinput()"); }
	}
}
/*-------------------------------------------------------------------------*/
void loudclickinput(int n)
/*-------------------------------------------------------------------------*/
{
switch(n) {
	case 0: { calibrating=TRUE; calibrate(10.0,n); break;} 
	case 1: { /* observation checked -- Figure 9 in Ruggero 1995 */
		/* animal L110, time-domain clicks */
        	data[n].signal=signaltype=CLICK;
        	data[n].dB=dB[0]=92;
		level[0]=scale*exp10(dB[0]/20.0);  
		observe=10.0;
		intone=0;
        	data[n].freq=0.0;
		data[n].nfeature=5;
		/* feature is time in ms -- onset time */ 
        	data[n].feature[0]=.2; 
		/* feature is time in ms -- decay */ 
        	data[n].feature[1]=1.3; 
		/* feature is time in ms -- duration of spindle */ 
        	data[n].feature[2]=1.6; 
		/* feature is peak-to-peak velocity in micrometers/sec */ 
        	data[n].feature[3]=3000.0; 
		/* feature is number of cycles to near-maximum */ 
        	data[n].feature[4]=1; 
		break;
		}
	case 2: { /* observation checked -- Figure 9 in Ruggero 1995 */
		/* animal L110, time-domain clicks */
        	data[n].signal=signaltype=CLICK;
        	data[n].dB=dB[0]=62;
		level[0]=scale*exp10(dB[0]/20.0);  
		observe=10.0;
		intone=0;
        	data[n].freq=0.0;
		data[n].nfeature=5;
		/* feature is time in ms -- onset time */ 
        	data[n].feature[0]=.5; 
		/* feature is time in ms -- decay */ 
        	data[n].feature[1]=.5; 
		/* feature is time in ms -- duration of spindle */ 
        	data[n].feature[2]=1.3; 
		/* feature is peak-to-peak velocity in micrometers/sec */ 
        	data[n].feature[3]=200.0; 
		/* feature is number of cycles to near-maximum */ 
        	data[n].feature[4]=4; 
		break;
		}
	default: { error("input case undefined","loudclickinput()"); }
	}
}
/*-------------------------------------------------------------------------*/
