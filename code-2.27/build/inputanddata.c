#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gold.h>

/*-------------------------------------------------------------------------*/
void defineinput(
                 int n, 
                 int inptype,
                 int datatype)
/*-------------------------------------------------------------------------*/
{
    int ifreq,idB;

    scale= 2.0*1e-07;	 	/* 2 x 1e-06 mg/mm*sec2 */

    switch(inptype) 
    {
	    case diepInput: case viergeverInput: 
            {
		    /* for matching data from fig 5.2.3 in viergever's thesis
		       the input signal is a pure tone with 1 KHz frequency
		       and SPL not given - assume 10 dB (doesn't matter, its
		       linear and axis for SPL becomes arbitrarily zeroed to
		       lowest response*/
		    signaltype = TONE;
		    dB[0] =10.0;
		    freq[0]=1.0;
		    intone=1;
		    level[0]=scale*exp10(dB[0]/20.0);  
		    nscanmodel=1;
		    break;
		    }
	    case userInput: 
            {
		    /* user sets the input */
		    signaltype = TONE;
		    idB = query("Input the decibels",1,100);
		    dB[0] =(double) idB;
		    ifreq = query("Input the frequency in Hz",100,100000);
		    freq[0]=(double) ifreq/1000;
		    intone=1;
		    level[0]=scale*exp10(dB[0]/20.0);  
		    nscanmodel=1;
		    break;
		    }
	    case fixedInput: 
            {
		    /* set the input to whatever */
		    signaltype = TONE;
		    dB[0] = 40.0;
		    freq[0]= 8.5;
		    intone=1;
		    level[0]=scale*exp10(dB[0]/20.0);  
		    nscanmodel=1;
		    break;
		    }
	    case neelyInput: 
            {
		    signaltype=TONE;
		    dB[0]=40.0;
		    freq[0]=1.60;
		    switch(n) 
            {
		        case 0: {
			        dB[0]=40.0;
			        break; }
		        case 1: {
			        dB[0]=.10;
			        break; }
		        case 2: {
			        dB[0]=100.0;
			        break; }
                default: {error("In inputanddata.c, defineinput()","neelyInput case undefined.");}
		    }
		    intone=1;
		    level[0]=scale*exp10(dB[0]/20.0);  
		    nscanmodel=1;
		    break;
		    }
	    case scanInput: 
            {
	        calibrating=FALSE;
	        data[n].inputnum=n;
	        data[n].dattype=datatype;
	        switch(datatype) 
            {
		        case liveclick: {nscanmodel=numlclick; break;}
		        case isointensity: {nscanmodel=numliso; break;}
		        case ttoned: {nscanmodel=numttoned; break;}
		        case ttones: {nscanmodel=numttones; break;}
		        case livetune: {nscanmodel=numltune; break;}
		        case quicklivetune: {nscanmodel=numqtune; break;}
		        case livegain: {nscanmodel=numgain; break;}
		        case loudclick: {nscanmodel=numloudclick; break;}
		        case loudiso: {nscanmodel=numloudiso; break;}
		        case postclick: {nscanmodel=numpclick; break;}
		        case posttune: {nscanmodel=numptune; break;}
		        case furoclick: {nscanmodel=numfclick; break;}
		        case furoinout: {nscanmodel=numfio; break;}
		        case furoclickfft: {nscanmodel=numfclickfft; break;}
		        default:{ error("datatype is undefined","in defineinput()"); }
		    }
	        switch(animal_status) 
            {
		        case live: {livescan(n,datatype); break; }
		        case furo: {furoscan(n,datatype); break; } 
		        case loud: {loudscan(n,datatype); break;}
		        case post: {postscan(n,datatype); break;}
		        default:{ error("Status is undefined","in defineinput()"); }
		    }	
	        data[n].observe=observe;
	        if ((!calibrating)&&(n!=0)) data[n].obsindex=data[n-1].obsindex;
	        if ((!calibrating)&&(n!=0)) data[n].fundindex=data[n-1].fundindex;
	        break;
	    } 
	    default:
        {
		    error("In defineinput.c", "Input not defined!");
		    break;
        }
	}
}
/*-------------------------------------------------------------------------*/
void calibrate(
               double obsfreq,
               int n)
/*-------------------------------------------------------------------------*/
{
	signaltype=TONE;
	dB[0]=40.0;
	freq[0]=obsfreq;
	observe=obsfreq;
	intone=1;
	level[0]=scale*exp10(dB[0]/20.0);  
	data[n].signal=signaltype;
    data[n].dB=dB[0];
    data[n].freq=freq[0];
    data[n].nfeature=0;
    data[n].cal=TRUE;
}
/*-------------------------------------------------------------------------*/
