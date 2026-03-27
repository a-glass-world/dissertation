#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"

/*-------------------------------------------------------------------------*/
int buildnet(
             int modeltype,
             int starting)
/*-------------------------------------------------------------------------*/
{
    int obj,chan,stability;

    stability=stable;
    buildmesh();
    buildgeometry();
    stability = buildimpedance();
    buildcoupling();

    /* do these builds in order, there are dependencies! */
    for (chan = 0;chan < nchannel;chan++) 
    {
	    if (channel[chan])  buildGamma(chan);
    }

    for (obj = 0;obj < nobject;obj++) 
    {
	    if ((fluidbound[obj]) && (object[obj]))
        {
            if (starting) initvgamma(obj);
		    else buildvgamma(obj);
        }
    }

    buildalpha(modeltype);
    buildA(FALSE);
    buildmisc();

    return(stability);
}
/*-------------------------------------------------------------------------*/
void buildmesh()
/*-------------------------------------------------------------------------*/
{
    int i;

    for (i=0;i<=length;i++) 
    {
	    mesh[i]=(double) i*cl/length;
	    if (i > 0) delta[i-1]=mesh[i]-mesh[i-1];
    }
} 
/*-------------------------------------------------------------------------*/
void buildgeometry() 
/*-------------------------------------------------------------------------*/
{
    int i,j;
    double h;

    if (gaindefined) 
        for (i = 0;i < length;i++) 
            gain[i] = gain_a11 * exp(gain_a12 *  mesh[i]);

    for (i=0;i<nobject;i++) 
    {
        if (object[i]) 
        {
            if ((plate[i]) || (fluidbound[i])) 
            { 
                for (j = 0;j < length;j++)  
                {
	                beta[i][j] = width[i][0] * exp(width[i][1] * mesh[j]);
                }
            }
        }
    }


    for (i = 0;i < nchannel;i++) 
    { 
        if (channel[i]) 
        {
	        for (j = 0;j < length;j++) 
            {
	            h = height[i][0] + height[i][1] * exp(height[i][2] * mesh[j] * mesh[j]);
	            area[i][j] = M_PI * (h/2.0) * (h/2.0);
            }
        }
    }
}
/*-------------------------------------------------------------------------*/
int buildimpedance() 
/*-------------------------------------------------------------------------*/
{
    int i,j;
    int stability;

    stability=stable;
    for (i=0;i<nobject;i++) 
    {
        if ((object[i]) && (i != ELMO)) 
        {
	        for (j=0;j<length;j++) 
            {
     	        R[i][j]=resist[i][0]*exp(resist[i][1]*mesh[j]);
     	        S[i][j]=stiff[i][0]*exp(stiff[i][1]*mesh[j]);
     	        M[i][j]=mass[i][0]*exp(mass[i][1]*mesh[j]);
		        damp[i][j]=R[i][j]/M[i][j];
		        lambda[i][j]=S[i][j]/M[i][j];
		        if (i==OHC) 
                {
			        maxld[j]=lambda[OHC][j];
  	               	        minld[j]=lambda[OHC][j]/10.0;
    		                rangeld[j]=maxld[j] - minld[j];
		        }
	        }
        }
    }

    lomap=fpx(S[BAM][length-1],M[BAM][length-1],R[BAM][length-1],&stability)/1000;
    himap=fpx(S[BAM][0],M[BAM][0],R[BAM][0],&stability)/1000;

    return(stability);
}
/*-------------------------------------------------------------------------*/
void buildcoupling() 
/*-------------------------------------------------------------------------*/
{
    int i,j,k;

    if (!CUTE) return;

    for (i=0;i<nobject;i++) 
    {
        for (j=0;j<nobject;j++) 
        {
	        if ((object[i])&&(object[j])&&(coupled[i][j]))
            {
	            for (k=0;k<length;k++) 
                {
		            cR[i][j][k] =  cresist[i][j][0]*exp(cresist[i][j][1]* mesh[k]);
		            cS[i][j][k] =  cstiff[i][j][0]*exp(cstiff[i][j][1]* mesh[k]);
		        }
            }
	    }
    }

    for (i=0;i<length;i++)
        staticS[i] =  cstiff[TM][ACRL][0]*exp(cstiff[TM][ACRL][1]* mesh[i]);

}
/*-------------------------------------------------------------------------*/
void buildalpha(int modeltype)
/*-------------------------------------------------------------------------*/
{
    int i,j,k;
    double maxalpha;

    maxalpha=50.0;

    for (i=0;i<nobject;i++) 
        for (j=0;j<length;j++) 
            alpha[i][j]=0.0;


    /* for checking against others, force it. */
    switch(modeltype) 
    {
        case neely2DOF: 
        {
		    for (k=0;k<length;k++) 
            {
			    alpha[BAM][k]=2.0/M[BAM][k];
			    if (fluidbound[TM]) alpha[TM][k]=2.0/M[TM][k];
			}
	    	break;
	    }
	    case diependaal1DOF: 
        {
		    for (k=0;k<length;k++) alpha[BAM][k]=.40;
		    break;
		}
	    default:
        {
		    if (MULTI) 
            { 
		        for (i=0;i<nobject;i++) 
                { 
			        if ((fluidbound[i])&&(object[i])) 
                    {
			            for (k=0;k<length;k++) 
                        {
			                alpha[i][k]=density*(beta[i][k]/M[i][k])
			                     *(1/(Gamma0[UC[i]][k]*area[UC[i]][k]) 
			                      +1/(Gamma0[LC[i]][k]*area[LC[i]][k])); 
                        }
                    }
                }
                /* just to temporarily override */
                for (k=0;k<length;k++) 
                {
			        alpha[BAM][k]=2.0/M[BAM][k];
			        if (fluidbound[TM]) alpha[TM][k]=2.0/M[TM][k];
			    }

            }
		    else 
            { 
		        for (i=0;i<nobject;i++) 
                {
			        if ((fluidbound[i])&&(object[i])) 
                    {
			            for (k=0;k<length;k++) 
                        {
			                alpha[i][k]=2.0*density*(beta[i][k]/M[i][k])
                                *(1/(Gamma0[SM][k]*area[SM][k]));
                        }
                    }
                }
            }
		    break;
		}
    }
} 
/*-------------------------------------------------------------------------*/
void buildmisc() 
/*-------------------------------------------------------------------------*/
{
    int x;
    double b0 = 1.45;

    ZweigN = sqrt(cl*cl*beta[BAM][0]/mass[BAM][0])/16.0;
    expb0 = exp(b0);

    for (x=0;x<length;x++) 
    {
                /* calculate basolateral resistance for OHC at location x*/
                meshRb[x] = Rb0 * exp(Rb1*mesh[x]); /* MOhms */
                /* calculate capacitance for OHC at location x*/
                meshCa0[x]=Ca0 * exp(Ca1*mesh[x]);
                meshCb0[x]=Cb0 * exp(Cb1*mesh[x]);
    }

}
