#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gold.h>


/*-------------------------------------------------------------------------*/
void livescan(
              int n,
              int measuretype)
/*-------------------------------------------------------------------------*/
{
static int started, hold;
int in;

if (!started) { hold=n; started=TRUE;}

in=n-hold;
switch(measuretype) {
	case liveclick: {liveclickinput(in); break;}
	case isointensity: {liveisoinput(in); break;}
	case ttoned: {livettonedinput(in);break;}
	case ttones: {livettonesinput(in);break;}
	case livetune: {tuningcurve(in); break;}
	case quicklivetune: {quicktune(in); break;}
	case livegain: {gaininput(in); break;}
	default: { error("measure type is not defined","livescan()");}
	}

/* 
for (i=0;i<data[in].nfeature;i++) data[in].weight[i]=1.0;
data[in].nweight=data[in].nfeature;
weights are all equal now... */

}
/*-------------------------------------------------------------------------*/
void liveclickinput(int n)
/*-------------------------------------------------------------------------*/
{

switch(n) {
	case 0: { calibrating=TRUE; calibrate(9.0,n); break;}
	case 1: { /* Live L13 from Nola - sequence 37 */
		data[n].signal=signaltype=CLICK;
        	data[n].dB=dB[0]=26;
		level[0]=scale*exp10(dB[0]/20.0);  
		observe=9.0;
		intone=0;
        	data[n].freq=0.0;
		data[n].nfeature=5;
		/* feature is time in ms -- onset time */ 
        	data[n].feature[0]=.8; 
		/* feature is time in ms -- decay */ 
        	data[n].feature[1]=1.3; 
		/* feature is time in ms -- duration of spindle */ 
        	data[n].feature[2]=0.4; 
		/* feature is peak-to-peak velocity in micrometers/sec */ 
        	data[n].feature[3]=145; 
		/* feature is number of cycles to near-maximum */ 
        	data[n].feature[4]=5; 
		break;
		}
	case 2: { /* Live L13 from Nola - sequence 37 */
		data[n].signal=signaltype=CLICK;
        	data[n].dB=dB[0]=56;
		level[0]=scale*exp10(dB[0]/20.0);  
		observe=9.0;
		intone=0;
        	data[n].freq=0.0;
		data[n].nfeature=5;
		/* feature is time in ms -- onset time */ 
        	data[n].feature[0]=.6; 
		/* feature is time in ms -- decay */ 
        	data[n].feature[1]=3.0; 
		/* feature is time in ms -- duration of spindle */ 
        	data[n].feature[2]=0.4; 
		/* feature is peak-to-peak velocity in micrometers/sec */ 
        	data[n].feature[3]=800; 
		/* feature is number of cycles to near-maximum */ 
        	data[n].feature[4]=4; 
		break;
		}
	case 3: { /* Live L13 from Nola - sequence 37 */
		data[n].signal=signaltype=CLICK;
        	data[n].dB=dB[0]=86;
		level[0]=scale*exp10(dB[0]/20.0);  
		observe=9.0;
		intone=0;
        	data[n].freq=0.0;
		data[n].nfeature=5;
		/* feature is time in ms -- onset time */ 
        	data[n].feature[0]=.3; 
		/* feature is time in ms -- decay */ 
        	data[n].feature[1]=5.4; 
		/* feature is time in ms -- duration of spindle */ 
        	data[n].feature[2]=0.4; 
		/* feature is peak-to-peak velocity in micrometers/sec */ 
        	data[n].feature[3]=2000; 
		/* feature is number of cycles to near-maximum */ 
        	data[n].feature[4]=3; 
		break;
		}
	default: { error("input case undefined","liveclickinput()"); }
	}
}
/*-------------------------------------------------------------------------*/
void liveisoinput(int n)
/*-------------------------------------------------------------------------*/
{

switch(n) {
	case 0: { calibrating=TRUE; calibrate(9.0,n); break;}
	case 1: { /* observation checked -- for low levels Figure 1.11 */
		    /* recalibrating for higher levels is cheating */
		    /* animal L13, Isointensity curve for 30 dB at 7 KHz */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=30;
		    level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=7.0;
		    observe=9.0;
		    intone=1;
		    data[n].nfeature=1;
		    /* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=15; 
		    break;
		    }
	case 2: { /* animal L13, Isointensity curve for 30 dB at 8 KHz */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=30;
		    level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=8.0;
		    observe=9.0;
		    intone=1;
		    data[n].nfeature=1;
		    /* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=215;
        	break;
		    }
	case 3: { /* animal 13, Isointensity curve for 30 dB at 9 KHz */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=30;
		    level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=9.0;
		    observe=9.0;
		    intone=1;
        	data[n].nfeature=1;
		    /* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=230;
        	break;
        	}
	case 4: { /* animal 13, Isointensity curve for 30 dB at 10 KHz */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=30;
		    level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=10.0;
		    observe=9.0;
		    intone=1;
        	data[n].nfeature=1;
		    /* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=80;
		    break;
		    }
	case 5: { /* animal 13, Isointensity curve for 50 dB at 7 KHz */
        	data[n].signal= signaltype=TONE;
        	data[n].dB= dB[0]=50;
		    level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=7.0;
		    observe=9.0;
		    intone=1;
        	data[n].nfeature=1;
		    /* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=180;
		    break;
		    }	
	case 6: { /* animal 13, Isointensity curve for 50 dB at 8 KHz */
        	data[n].signal= signaltype=TONE;
        	data[n].dB= dB[0]=50;
		    level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=8.0;
		    observe=9.0;
		    intone=1;
        	data[n].nfeature=1;
		    /* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=500;
		    break;
		    }	
	case 7: { /* animal 13, Isointensity curve for 50 dB at 9 KHz */
        	data[n].signal= signaltype=TONE;
        	data[n].dB= dB[0]=50;
		    level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=9.0;
		    observe=9.0;
		    intone=1;
        	data[n].nfeature=1;
		    /* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=400;
		    break;
	    	}
	case 8: { /* animal 13, Isointensity curve for 50 dB at 10 KHz */
        	data[n].signal= signaltype=TONE;
        	data[n].dB= dB[0]=50;
		    level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=10.0;
		    observe=9.0;
		    intone=1;
        	data[n].nfeature=1;
		    /* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=130;
		    break;
		    }
	case 9: { /* animal 13, Isointensity curve for 80 dB at 3 KHz */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=80;
	    	level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=3.0;
	    	observe=9.0;
		    intone=1;
        	data[n].nfeature=1;
		    /* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=50;
		    break;
		    }
	case 10: { /* animal 13, Isointensity curve for 80 dB at 7 KHz */
       	 	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=80;
	    	level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=7.0;
    		observe=9.0;
    		intone=1;
      		data[n].nfeature=1;
    		/* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=1500;
    		break;
    		}
	case 11: { /* animal 13, Isointensity curve for 80 dB at 8 KHz */
       	 	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=80;
	    	level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=8.0;
	    	observe=9.0;
	    	intone=1;
      		data[n].nfeature=1;
	    	/* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=1000;
	    	break;
        	}
	case 12: { /* animal 13, Isointensity curve for 80 dB at 9 KHz */
       	 	data[n].signal= signaltype=TONE;
        	data[n].dB=dB[0]=80;
		    level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=9.0;
		    observe=9.0;
		    intone=1;
      		data[n].nfeature=1;
		    /* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=600;
	    	break;
	        }
	case 13: { /* animal 13, Isointensity curve for 80 dB at 10 KHz */
       	 	data[n].signal= signaltype=TONE;
        	data[n].dB=dB[0]=80;
		    level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=10.0;
		    observe=9.0;
		    intone=1;
      		data[n].nfeature=1;
		    /* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=200;
		    break;
        	}
	default: {  error("input case undefined","liveisoinput()"); }

	}
}
/*-------------------------------------------------------------------------*/
void livettonedinput(int n)
/*-------------------------------------------------------------------------*/
{
double tempdB;
int i;

switch(n) {
/* in two tones make sure the first frequencies are the higher */
	case 0: { calibrating=TRUE; calibrate(8.5,n); break;}
	case 1: { /* checked observation. animal L17, 90 dB */ 
        	data[n].signal=signaltype=TTONED;
        	data[n].dB= dB[0]=90;
		    dB[1]=90;
        	data[n].freq=freq[0]=7.79;
		    freq[1]=7.08;
        	/* there are 6 peaks at
			 3f1 - f2, 2f1 - f2, f1, f2, 2f2 - f1(CF), 3f2 - 2f1 */
        	/* these features are measured in micrometers/sec */
        	data[n].nfeature=6;
        	data[n].feature[0]=25; /* 3f1 - 2f2 */
        	data[n].feature[1]=50; /* 2f1 - f2 */
        	data[n].feature[2]=3000; /* f1 */
        	data[n].feature[3]=1000; /* f2 */
        	data[n].feature[4]=100; /* 2f2 - f1 */
        	data[n].feature[5]=10; /* 3f2 - 2f1 */
        	for (i=0;i<data[n].nfeature;i++) {
                	tempdB=20.0*log10(data[n].feature[i]);
                	data[n].feature[i]=tempdB; }
		    observe=8.5;
		    intone=2;
		    level[0]=scale*exp10(dB[0]/20.0);  
		    level[1]=scale*exp10(dB[1]/20.0);  
		    break;
		    }
	case 2: { /* checked observation. animal L17, 70 dB */ 
        	data[n].signal=signaltype=TTONED;
        	data[n].dB= dB[0]=70;
		    dB[1]=70;
        	data[n].freq=freq[0]=7.79;
		    freq[1]=7.08;
        	/* there are 6 peaks at
			 3f1 - f2, 2f1 - f2, f1, f2, 2f2 - f1(CF), 3f2 - 2f1 */
        	/* these features are measured in micrometers/sec */
        	data[n].nfeature=6;
        	data[n].feature[0]=5; /* 3f1 - 2f2 */
        	data[n].feature[1]=100; /* 2f1 - f2 */
        	data[n].feature[2]=1000; /* f1 */
       	 	data[n].feature[3]=800; /* f2 */
      	  	data[n].feature[4]=110; /* 2f2 - f1 */
       		data[n].feature[5]=30; /* 3f2 - 2f1 */
        	for (i=0;i<data[n].nfeature;i++) {
                	tempdB=20.0*log10(data[n].feature[i]);
                	data[n].feature[i]=tempdB; }
		    observe=8.5;
		    intone=2;
		    level[0]=scale*exp10(dB[0]/20.0);  
		    level[1]=scale*exp10(dB[1]/20.0);  
		    break;
		    }
	case 3: { calibrating=TRUE; calibrate(8.0,n); break;}
	case 4: { /* checked observation. animal L26, 80 dB */ 
        	data[n].signal=signaltype=TTONED;
        	data[n].dB= dB[0]=80;
        	data[n].freq=freq[0]=10.8;
        	/* there are 6 peaks at
			 3f1 - f2, 2f1 - f2, f1, f2, 2f2 - f1(CF), 3f2 - 2f1 */
        	/* these features are measured in micrometers/sec */
        	data[n].nfeature=6;
        	data[n].feature[0]=30; /* 3f1 - 2f2 */
        	data[n].feature[1]=40; /* 2f1 - f2 */
        	data[n].feature[2]=35; /* f1 */
        	data[n].feature[3]=38; /* f2 */
        	data[n].feature[4]=0; /* 2f2 - f1 */
        	data[n].feature[5]=0; /* 3f2 - 2f1 */
		    dB[1]=80;
		    freq[1]=9.41;
		    observe=8.0;
		    intone=2;
		    level[0]=scale*exp10(dB[0]/20.0);  
		    level[1]=scale*exp10(dB[1]/20.0);  
		    break;
		    }
	case 5: { /* checked observation. animal L26, 40 dB */ 
        	data[n].signal=signaltype=TTONED;
        	data[n].dB=dB[0]=40;
        	data[n].freq= freq[0]=10.8;
        	/* there are 6 peaks at
			 3f1 - f2, 2f1 - f2, f1, f2, 2f2 - f1(CF), 3f2 - 2f1 */
        	/* these features are measured in micrometers/sec */
        	data[n].nfeature=6;
        	data[n].feature[0]=0; /* 3f1 - 2f2 */
        	data[n].feature[1]=23; /* 2f1 - f2 */
        	data[n].feature[2]=20; /* f1 */
        	data[n].feature[3]=0; /* f2 */
        	data[n].feature[4]=0; /* 2f2 - f1 */
        	data[n].feature[5]=0; /* 3f2 - 2f1 */
		    dB[1]=40;
		    freq[1]=9.41;
		    observe=8.0;
		    intone=2;
		    level[0]=scale*exp10(dB[0]/20.0);  
		    level[1]=scale*exp10(dB[1]/20.0);  
		    break;
		    }
	default: { error("input case undefined","livettonedinput()"); }
	}
}
/*-------------------------------------------------------------------------*/
void livettonesinput(int n)
/*-------------------------------------------------------------------------*/
{

switch(n) {
	case 0: { calibrating=TRUE; calibrate(8.0,n); break;}
	case 1: { /* animal L29, suppressor 10 KHz, Probe 8 KHz */ 
		observe=8.0;    /* observe at the probe */
		signaltype=TTONES;
		dB[0]=70;
		dB[1]=40;
		freq[0]=10.0;   /* suppressor */
		freq[1]=8.0; 	/* probe */
		intone=2;
		level[0]=scale*exp10(dB[0]/20.0);  
		level[1]=scale*exp10(dB[1]/20.0);  
        	data[n].signal=signaltype;
		data[n].freq=freq[0];
		data[n].dB=dB[0];
        	data[n].nfeature=1;
		/* velocity in micrometers per sec */
		data[n].feature[0]=15; 
		break;
		}
	case 2: {
		signaltype=TTONES;
		dB[0]=70;
		dB[1]=70;
		freq[0]=10.0;
		freq[1]=8.0;
		observe=8.0; /* observe at the probe */
		intone=2;
		level[0]=scale*exp10(dB[0]/20.0);  
		level[1]=scale*exp10(dB[1]/20.0);  
        	data[n].signal=signaltype;
		data[n].freq=freq[0];
		data[n].dB=dB[0];
        	data[n].nfeature=1;
		/* velocity in micrometers per sec */
		data[n].feature[0]=500; 
		break;
		}
	case 3: { /* this is pure tone 8.0 KHz for comparison*/
		signaltype=TONE;
		freq[0]=8.0;
		dB[0]=40;
		observe=8.0; /* observe at the probe */
		intone=1;
		level[0]=scale*exp10(dB[0]/20.0);  
		data[n].signal=signaltype;
		data[n].freq=freq[0];
		data[n].dB=dB[0];
        	data[n].nfeature=1;
		/* velocity in micrometers per sec */
		data[n].feature[0]=70; 
		break;
		}
	case 4: { /* this is pure tone at probe 8.0 KHz */
		signaltype=TONE;
		dB[0]=70;
		freq[0]=8.0;
		observe=8.0; /* observe at the probe */
		intone=1;
		level[0]=scale*exp10(dB[0]/20.0);  
		data[n].signal=signaltype;
		data[n].freq=freq[0];
		data[n].dB=dB[0];
        	data[n].nfeature=1;
		/* velocity in micrometers per sec */
		data[n].feature[0]=600; 
		break;
		}
	case 5: { calibrating=TRUE; calibrate(6.8,n); break;}

	case 6: { /* Suppressor at 1 KHz, 85 dB, Probe at 6.8 KHz */
		signaltype=TTONES;
		dB[0]=70;
		dB[1]=85;
		freq[0]=6.8; /* probe */
		freq[1]=1.0; /* low side suppressor */
		observe=6.8;
		intone=2;
		level[0]=scale*exp10(dB[0]/20.0);  
		level[1]=scale*exp10(dB[1]/20.0);  
		data[n].signal=signaltype;
		data[n].freq=freq[0];
		data[n].dB=dB[0];
        	data[n].nfeature=1;
		/* velocity in micrometers per sec */
		data[n].feature[0]=300; 
		break;
		}
	case 7: { /* Suppressor at 1 KHz, 85 dB, Probe at 6.8 KHz */
		observe=6.8;
		signaltype=TTONES;
		dB[0]=85;
		dB[1]=85;
		freq[0]=6.8;
		freq[1]=1.0;
		intone=2;
		level[0]=scale*exp10(dB[0]/20.0);  
		level[1]=scale*exp10(dB[1]/20.0);  
		data[n].signal=signaltype;
		data[n].freq=freq[0];
		data[n].dB=dB[0];
        	data[n].nfeature=1;
		/* velocity in micrometers per sec */
		data[n].feature[0]=900; 
		break;
		}
	case 8: { /* this is pure tone at probe 6.8 KHz */
		signaltype=TONE;
		dB[0]=70;
		freq[0]=6.8;
		observe=6.8; /* observe at the probe */
		intone=1;
		level[0]=scale*exp10(dB[0]/20.0);  
		data[n].signal=signaltype;
		data[n].freq=freq[0];
		data[n].dB=dB[0];
        	data[n].nfeature=1;
		/* velocity in micrometers per sec */
		data[n].feature[0]=700; 
		break;
		}
	case 9: { /* this is pure tone at probe 6.8 KHz */
		signaltype=TONE;
		dB[0]=85;
		freq[0]=6.8;
		observe=6.8; /* observe at the probe */
		intone=1;
		level[0]=scale*exp10(dB[0]/20.0);  
		data[n].signal=signaltype;
		data[n].freq=freq[0];
		data[n].dB=dB[0];
        	data[n].nfeature=1;
		/* velocity in micrometers per sec */
		data[n].feature[0]=1200; 
		break;
		}	
	default: { error("input case undefined","livettonesinput()"); }
	}	
}/*-------------------------------------------------------------------------*/
void tuningcurve(int n)
/*-------------------------------------------------------------------------*/
{

switch(n) {
	case 0: 
        { 
            calibrating=TRUE; calibrate(8.0,n); 
            break;
        }
	case 1:
        {
            /* observation checked -- Figures 6 and 12 in Robles 1986 */
		    /* animal M044, Isoresponse curves */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=90;
		    level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=1.0;
		    observe=8.0;
		    intone=1;
		    data[n].nfeature=2;
		    /* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=100.0; 
        	data[n].feature[1]=-1.4;
		    break;
		}
	case 2:
        {
            /* observation checked -- Figures 6 and 12 in Robles 1986 */
		    /* animal M044, Isoresponse curves */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=88;
		    level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=2.0;
		    observe=8.0;
		    intone=1;
		    data[n].nfeature=2;
		    /* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=100.0; 
		    /* feature is phase */ 
        	data[n].feature[1]=-1.5;
		    break;
		}
	case 3:
        {
            /* observation checked -- Figures 6 and 12 in Robles 1986 */
		    /* animal M044, Isoresponse curves */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=79;
		    level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=3;
		    observe=8.0;
		    intone=1;
		    data[n].nfeature=2;
		    /* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=100.0; 
		    /* feature is phase in micrometers/sec */ 
        	data[n].feature[1]=-1.0;
		    break;
		    }	
	case 4:{/* observation checked -- Figures 6 and 12 in Robles 1986 */
		/* animal M044, Isoresponse curves */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=68;
	    	level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=4;
	    	observe=8.0;
	    	intone=1;
	    	data[n].nfeature=2;
	    	/* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=100.0; 
	    	/* feature is phase in micrometers/sec */ 
        	data[n].feature[1]=-1.4;
	    	break;
	    	}	
	case 5:{/* observation checked -- Figures 6 and 12 in Robles 1986 */
		    /* animal M044, Isoresponse curves */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=56;
		    level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=5;
		    observe=8.0;
		    intone=1;
		    data[n].nfeature=2;
		    /* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=100.0; 
		    /* feature is phase in micrometers/sec */ 
        	data[n].feature[1]=-2.0;
		    break;
		    }	
	case 6:{/* observation checked -- Figures 6 and 12 in Robles 1986 */
		/* animal M044, Isoresponse curves */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=48;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=6;
		observe=8.0;
		intone=1;
		data[n].nfeature=2;
		/* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=100.0; 
		/* feature is phase in micrometers/sec */ 
        	data[n].feature[1]=-2.3;
		break;
		}	
	case 7:
        {
            /* observation checked -- Figures 6 and 12 in Robles 1986 */
		    /* animal M044, Isoresponse curves */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=30;
		    level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=7;
		    observe=8.0;
		    intone=1;
		    data[n].nfeature=2;
		    /* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=100.0; 
		    /* feature is phase in micrometers/sec */ 
        	data[n].feature[1]=-2.8;
		    break;
		}	
	case 8:{/* observation checked -- Figures 6 and 12 in Robles 1986 */
		/* animal M044, Isoresponse curves */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=15;
		    level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=8;
		    observe=8.0;
		    intone=1;
		    data[n].nfeature=2;
		    /* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=100.0; 
		    /* feature is phase in micrometers/sec */ 
        	data[n].feature[1]=-4.3;
		    break;
		    }	
	case 9:{/* observation checked -- Figures 6 and 12 in Robles 1986 */
		/* animal M044, Isoresponse curves */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=12;
		    level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=8.5;
		    observe=8.0;
		    intone=1;
		    data[n].nfeature=2;
		    /* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=100.0; 
		    /* feature is phase in micrometers/sec */ 
        	data[n].feature[1]=-5.2;
		    break;
		    }	
	case 10:{/*observation checked -- Figures 6 and 12 in Robles 1986 */
		    /* animal M044, Isoresponse curves */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=18;
		    level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=9;
		    observe=8.0;
		    intone=1;
		    data[n].nfeature=2;
		    /* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=100.0; 
		    /* feature is phase in micrometers/sec */ 
        	data[n].feature[1]=-6.2;
		break;
		}	
	case 11:{/*observation checked -- Figures 6 and 12 in Robles 1986 */
		/* animal M044, Isoresponse curves */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=40;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=9.5;
		observe=8.0;
		intone=1;
		data[n].nfeature=2;
		/* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=100.0; 
		/* feature is phase in micrometers/sec */ 
        	data[n].feature[1]=-7.0;
		break;
		}	
	case 12:{/*observation checked -- Figures 6 and 12 in Robles 1986 */
		/* animal M044, Isoresponse curves */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=60;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=10;
		observe=8.0;
		intone=1;
		data[n].nfeature=2;
		/* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=100.0; 
		/* feature is phase in micrometers/sec */ 
        	data[n].feature[1]=-8.0;
		break;
		}	
	case 13:{/*observation checked -- Figures 6 and 12 in Robles 1986 */
		/* animal M044, Isoresponse curves */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=76;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=10.5;
		observe=8.0;
		intone=1;
		data[n].nfeature=2;
		/* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=100.0; 
		/* feature is phase in micrometers/sec */ 
        	data[n].feature[1]=-9.2;
		break;
		}	
	case 14:{/*observation checked -- Figures 6 and 12 in Robles 1986 */
		/* animal M044, Isoresponse curves */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=82;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=11;
		observe=8.0;
		intone=1;
		data[n].nfeature=2;
		/* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=100.0; 
		/* feature is phase in micrometers/sec */ 
        	data[n].feature[1]=-9.0; 
		break;
		}	
	default: { error("input case undefined","tuningcurve()"); }
	}
}
/*-------------------------------------------------------------------------*/
void quicktune(int n)
/*-------------------------------------------------------------------------*/
{

switch(n) {
	case 0: 
        { 
            calibrating=TRUE; calibrate(8.5,n); 
            break;
        }
	case 1:
        {
            /* observation checked -- Figures 6 and 12 in Robles 1986 */
		    /* animal M044, Isoresponse curves */
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=90;
		    level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=1.0;
		    observe=8.0;
		    intone=1;
		    data[n].nfeature=2;
		    /* feature is velocity in micrometers/sec */ 
        	data[n].feature[0]=500.0; 
        	data[n].feature[1]=-1.4;
		    break;
		}
	case 2:{
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=56;
		    level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=5;
		    observe=8.0;
		    intone=1;
		    data[n].nfeature=2;
		   
        	data[n].feature[0]=100.0; 
		 
        	data[n].feature[1]=-2.0;
		    break;
		    }	
	case 3:{
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=12;
		    level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=8.5;
		    observe=8.0;
		    intone=1;
		    data[n].nfeature=2;
		
        	data[n].feature[0]=200.0; 
		  
        	data[n].feature[1]=-5.2;
		    break;
		    }	
	case 4:{
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=82;
		    level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=11;
		    observe=8.0;
		    intone=1;
		    data[n].nfeature=2;

        	data[n].feature[0]=100.0; 
	
        	data[n].feature[1]=-9.0; 
		    break;
		    }
       
	default: { error("input case undefined","tuningcurve()"); }
	}
}
/*-------------------------------------------------------------------------*/
void gaininput(int n)
/*-------------------------------------------------------------------------*/
{
switch(n) {
/* animal L113, Gain for 12,14.5, 21, 30,40,50,60,70,80,90,100,108 dB */
/* observations is at CP for 10 KHz, input is always 10 KHz */
/* feature is gain = velocity (micrometers/sec)/(.02*10^(dB/20) */ 
	case 0: { calibrating=TRUE; calibrate(10.0,n); break;}
	case 1:{
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=12;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=10.0;
		observe=10.0;
		intone=1;
		data[n].nfeature=1;
        	data[n].feature[0]=858.46; 
		break; }
	case 2:{
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=14.5;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=10.0;
		observe=10.0;
		intone=1;
		data[n].nfeature=1;
        	data[n].feature[0]=912.76;
		break; }
	case 3:{
        	data[n].signal=signaltype=TONE;
        	data[n].dB=dB[0]=21.0;
		    level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=10.0;
		    observe=10.0;
		    intone=1;
		    data[n].nfeature=1;
        	data[n].feature[0]=793.48;
		    break; }
	case 4:{
        	data[n].dB=dB[0]=30.0;
        	data[n].feature[0]=470.61;
        	data[n].signal=signaltype=TONE;
		    level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=10.0;
		    observe=10.0;
		    intone=1;
		    data[n].nfeature=1;
		    break; }
	case 5:{
        	data[n].dB=dB[0]=40.0;
        	data[n].feature[0]=323.10;
        	data[n].signal=signaltype=TONE;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=10.0;
		observe=10.0;
		intone=1;
		data[n].nfeature=1;
		break; }
	case 6:{
        	data[n].dB=dB[0]=50.0;
        	data[n].feature[0]=102.70;
        	data[n].signal=signaltype=TONE;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=10.0;
		observe=10.0;
		intone=1;
		data[n].nfeature=1;
		break; }
	case 7:{
        	data[n].dB=dB[0]=60.0;
        	data[n].feature[0]=43.13;
        	data[n].signal=signaltype=TONE;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=10.0;
		observe=10.0;
		intone=1;
		data[n].nfeature=1;
		break; }
	case 8:{
        	data[n].dB=dB[0]=70.0;
        	data[n].feature[0]=18.25;
        	data[n].signal=signaltype=TONE;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=10.0;
		observe=10.0;
		intone=1;
		data[n].nfeature=1;
		break; }
	case 9:{
        	data[n].dB=dB[0]=80.0;
        	data[n].feature[0]=7.48;
        	data[n].signal=signaltype=TONE;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=10.0;
		observe=10.0;
		intone=1;
		data[n].nfeature=1;
		break; }
	case 10:{
        	data[n].dB=dB[0]=90.0;
        	data[n].feature[0]=2.41;
        	data[n].signal=signaltype=TONE;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=10.0;
		observe=10.0;
		intone=1;
		data[n].nfeature=1;
		break; }
	case 11:{
        	data[n].dB=dB[0]=100.0;
        	data[n].feature[0]=.58;
        	data[n].signal=signaltype=TONE;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=10.0;
		observe=10.0;
		intone=1;
		data[n].nfeature=1;
		break; }
	case 12:{
        	data[n].dB=dB[0]=108.0;
        	data[n].feature[0]=.31;
        	data[n].signal=signaltype=TONE;
		level[0]=scale*exp10(dB[0]/20.0);  
        	data[n].freq=freq[0]=10.0;
		observe=10.0;
		intone=1;
		data[n].nfeature=1;
		break; }
	default: { error("input case undefined","gaininput()"); }
	}	
}
/*-------------------------------------------------------------------------*/
