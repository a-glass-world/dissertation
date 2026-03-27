/*-------------------------------------------------------------------------*/
/* This file runs the fluid solution for the time-domain in 1D.		   */
/*-------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"

#define DEBUG2 FALSE 
#define OVERRIDE FALSE

/*-------------------------------------------------------------------------*/
int fluids(
           int i,
           double time,
           int in,
           int out,
           int DEBUG)
/*-------------------------------------------------------------------------*/
{
    int stability;

    /* fluid globals */
    YF=in;
    DF=out;
    iO=i;

    CalculateG();
    CalculateGc();
    if (fluidbound[iO]!=TRUE)
    {
        QuickDerive(); 
    }
    else 
    {
	    GetSignal(time);
	    totsignal=BoundaryConditions();
	    CalculateK(totsignal);
	    CalculateP();
	    Derive();
	    if (MIDEAR) 
            DeriveStapes(); 
	}

    stability=CheckUnstable(); 

    if (DEBUG) fluiddebug(time,YF,DF,iO); 

    return(stability);
}

/*-------------------------------------------------------------------------*/
void GetSignal(double time)
/*-------------------------------------------------------------------------*/
{
    double dBmag, A, Omega, Sinwt,x;
    int i;

    switch(signaltype) 
    {
        case PIP: case TONE: case TTONES: case TTONED: 
        {
		    finput=0.0;
		    for (i=0;i<intone;i++) 
            {
			    A=level[i];  
			    Omega=(double) 2.0*M_PI*time*freq[i];
			    Sinwt = sin(Omega);
			    finput+=A*Sinwt;
		    }
		    break;
	    }
        case CLICKFFT: case CLICK: 
        {
		    finput=0.0;
            /* both of these inputs have been tuned to give the correct */
            /* simulated response for the passive model  -- however, the nonRECIO */
            /* response is "better" - i.e. more linear */
		    if (RECIO) 
		    {
			    /* this is from the stapes data for ME0001 at 0 dBA */
			    x=getclick(time);
			    /* therefore do NOT use middleear ! */
			    if (MIDEAR) 
		    	    error("The click data is for the stapes",
		    	    "no middle ear transfer function should be used!");
			    /* x is dB stapes response for 105 dB (0 dBA) */
			    /* it scales linearly - Reccio, 1996 */
			    dBmag=(x/165)*dB[0];
			    /* the level is now slightly changed */
			    level[0]=5.0*scale*exp10(dBmag/20.0);
			    finput=level[0];
            }
		    else 
		    {
			    if (time<=clickdur) 
                    finput=120*level[0];
			    else 
                    level[0]=0.0;
            }
		    break;
        }
        default: 
	    {
		    error("in fluid.c, GetSignal():\n","Input signal undefined.."); 
		    break;
	    }
    }
}
/*-------------------------------------------------------------------------*/
void CalculateG()
/*-------------------------------------------------------------------------*/
{
    int i;

    for (i=0; i<length; i++) 
        G[i]= R[iO][i]*Y[YF][1][iO].loc[i] + S[iO][i]*Y[YF][0][iO].loc[i];

}
/*-------------------------------------------------------------------------*/
void CalculateGc() 
/*-------------------------------------------------------------------------*/
{
    int i, comp;

    if (!CUTE) return;

    /* if the stepsize is sufficiently small, compare couple using YF  */
    comp = YF;
    for (i=0;i<length;i++) 
    {
	    if ((coupledOHC) && (motileOHC)) 
            Gc[i] = getactivecouple(i,YF);
        else
            Gc[i] = getpassivecouple(i,YF);
	} 
}
/*------------------------------------------------------------------------*/
double getpassivecouple(
                 int i,
                 int comp)
/*------------------------------------------------------------------------*/
{
    int j;
    double x0, v0, x1, v1;
    double Rc = 0.0;
    double Sc = 0.0;

    x0 = Y[comp][0][iO].loc[i];
    v0 = Y[comp][1][iO].loc[i];

    for (j = 0;j < nobject;j++) 
    {
        if (coupled[iO][j]) 
        {
            if ((iO == j) && (i != 0))
            {
                x1=Y[comp][0][j].loc[i-1];
		        v1=Y[comp][1][j].loc[i-1];

            }
            else
            {
                x1=Y[comp][0][j].loc[i];
		        v1=Y[comp][1][j].loc[i];
            }

		    Rc+=(cR[iO][j][i])*(v0 - v1);
            Sc+=(cS[iO][j][i])*(x0 - x1);
     	}
    }

    return(Rc + Sc);

}

/*------------------------------------------------------------------------*/
double getactivecouple(
                 int i,
                 int comp)
/*------------------------------------------------------------------------*/
{
    int j;
    double Sc = 0.0;
    double Rc = 0.0;
    double x0,v0,x1,v1;
    double xlc,vlc;

    x0 = Y[comp][0][iO].loc[i];
    v0 = Y[comp][1][iO].loc[i];

    /* set up active couple */
    xlc = activeC[comp][0][OHC].loc[i]; 
    vlc = activeC[comp][1][OHC].loc[i]; 

    if (iO == OHC) 
    { 
	    /* only three couples possible .. */
	    if (coupled[OHC][BAM]) 
        {
		    /* in phase */
		    x1=Y[comp][0][BAM].loc[i];
		    v1=Y[comp][1][BAM].loc[i];
		    Rc+=(cR[iO][BAM][i]) * ((v0 + vlc) - v1);
		    Sc+=(cS[iO][BAM][i]) * ((x0 + xlc) - x1);
            if (OVERRIDE) 
            {
                /* override */
		        Rc+=(cR[iO][BAM][i]) * ((v0 + 2*vlc) - v1);
		        Sc+=(cS[iO][BAM][i]) * ((x0 + 2*xlc) - x1);
            }
		}

	    if (coupled[OHC][ACRL]) 
        {
		    /* anti-phase */
		    x1=Y[comp][0][ACRL].loc[i];
		    v1=Y[comp][1][ACRL].loc[i];
		    Rc+=(cR[iO][ACRL][i]) * ((v0 - vlc) - v1);
		    Sc+=(cS[iO][ACRL][i]) * ((x0 - xlc) - x1);
            if (OVERRIDE) 
            {
                /* override */
		        Rc+=(cR[iO][ACRL][i]) * ((v0) - v1);
		        Sc+=(cS[iO][ACRL][i]) * ((x0) - x1);
            }
		}
     
	    if (coupled[OHC][OHC]) 
        {
            /* combined motions cancel */
            if (i != 0) 
            {
		        x1 = Y[comp][0][OHC].loc[i-1];
		        v1 = Y[comp][1][OHC].loc[i-1];
		        Rc += (cR[iO][OHC][i]) * (v0 - v1);
		        Sc += (cS[iO][OHC][i]) * (x0 - x1);
            }
		} 
    }
    else 
    { 	
        for (j = 0; j < nobject; j++) 
        {
	        if (coupled[iO][j]) 
            {
		        x1 = Y[comp][0][j].loc[i];
		        v1 = Y[comp][1][j].loc[i];

		        if ((j == OHC) && (iO == ACRL)) 
                { 
			        /* anti-phase */
			        Rc += (cR[iO][j][i])*(v0 - (v1 - vlc));
			        Sc += (cS[iO][j][i])*(x0 - (x1 - xlc));
                    if (OVERRIDE) 
                    {
                        /* override */
			            Rc += (cR[iO][j][i])*(v0 - (v1));
			            Sc += (cS[iO][j][i])*(x0 - (x1));
                    }
		        }
		        else if ((j == OHC) && (iO == BAM)) 
                { 
			        /* in phase */
			        Rc += (cR[iO][j][i])*(v0 - (v1 + vlc));
			        Sc += (cS[iO][j][i])*(x0 - (x1 + xlc));
                    if (OVERRIDE) 
                    {
                        /* override */
			            Rc += (cR[iO][j][i])*(v0 - (v1 + 2*vlc));
			            Sc += (cS[iO][j][i])*(x0 - (x1 + 2*xlc));
                    }
		        }
		        else if ((j == ACRL) && (iO == ACRL))
                { 	
                    if (i != 0) {

                        x1 = Y[comp][0][j].loc[i-1];
		                v1 = Y[comp][1][j].loc[i-1];

			            Rc += (cR[iO][j][i])*(v0 - v1);
			            Sc += (cS[iO][j][i])*(x0 - x1);
                    }
     	        }
		        else 
                { 	
			        Rc += (cR[iO][j][i])*(v0 - v1);
			        Sc += (cS[iO][j][i])*(x0 - x1);
     	        }
	        }
        }
    }

    return(Rc + Sc);

}
/*-------------------------------------------------------------------------*/
double BoundaryConditions() 
/*-------------------------------------------------------------------------*/
{
    double fsignal, insignal;


    if (MIDEAR)
    {
	    Gme = T2m*Rm*Y[YF][1][iO].stapes + Sm*Y[YF][0][iO].stapes;
	    fsignal = density*(Gme - Tm*Am*finput)/Tm2mm_ms; 
    }
    else    
    {
        fsignal = finput; 
    }

    if (iO==BAM) 
        insignal=2.0*fsignal;
    else 
        insignal=0;

    return(insignal);
}
/*-------------------------------------------------------------------------*/
void QuickDerive()
/*-------------------------------------------------------------------------*/
{
    int i;

    for (i = 0;i < length;i++) 
    {
	    D[DF][0][iO].loc[i] = Y[YF][1][iO].loc[i];
	    D[DF][1][iO].loc[i] = (-Gc[i] - G[i])/M[iO][i];
    }
}
/*-------------------------------------------------------------------------*/
void CalculateK(double insignal)
/*-------------------------------------------------------------------------*/
{
    int i;

    /* the first term of K includes the boundary condition at the base */
    K[0] = 0.5 * alpha[iO][0] * (G[0] + Gc[0]) * delta[0] + insignal;
    for (i = 1; i < length; i++) 
	    K[i]=0.5 * alpha[iO][i] * (G[i] + Gc[i]) * (mesh[i+1] - mesh[i-1]);
}
/*-------------------------------------------------------------------------*/
void CalculateP()
/*-------------------------------------------------------------------------*/
{
    double dfY[length];
    int i;
    /*-----------------------------------------------------------------*/
    /* using lower diagonal L and upper diagonal U, solve  Y from      */
    /*		L*Y=K						                               */
    /*  and P from							                           */
    /*		U*P = Y	 		    		                               */
    /*-----------------------------------------------------------------*/

    dfY[0]=K[0]/L[iO][0][LDIAG];
    for (i=1;i<length;i++) 
	    dfY[i]=(K[i] - (L[iO][i][LSUBDIAG] * dfY[i-1])) / L[iO][i][LDIAG];

    P[DF].loc[length-1] = dfY[length-1]/U[iO][length-1][UDIAG];
    for (i = (length-2);i >= 0;i--) 
	    P[DF].loc[i] = 
            (dfY[i] - (U[iO][i][USUPRADIAG] * P[DF].loc[i+1])) / U[iO][i][UDIAG];
}
/*-------------------------------------------------------------------------*/
void DeriveStapes()
/*-------------------------------------------------------------------------*/
{
    /* this only works with fluid pressure in the scala media */
    D[DF][1][iO].stapes = (Tm*Am*finput+area[SM][0]*P[DF].loc[0]-Gme)/Tm2mm_ms;
    D[DF][0][iO].stapes = Y[YF][1][iO].stapes; 
}

/*-------------------------------------------------------------------------*/
void Derive()
/*-------------------------------------------------------------------------*/
{
int i;

for (i=0;i<length;i++) {
	D[DF][0][iO].loc[i]= Y[YF][1][iO].loc[i];
	D[DF][1][iO].loc[i]=(P[DF].loc[i]- G[i] - Gc[i])/M[iO][i];
	}
} 
/*-------------------------------------------------------------------------*/
int CheckUnstable()
/*-------------------------------------------------------------------------*/
{
    int i;
    double err;

    i=0;
    err=1.0;
    if ((Y[YF][0][iO].loc[i]>err)||(Y[YF][0][iO].loc[i]<(-err))) 
	{
        nameit(iO,stderr);
	    fprintf(stderr,"has unstable values %e, %e \n", 
		    Y[YF][0][iO].loc[i],Y[YF][0][iO].loc[i]);
	    
        if ((!estimating)&&(!scanning)&&(!ranging)) 
		    error("in fluids.c, CheckUnstable()","Unstable solution\n ");
	    else return(unstable);
	}

    return(stable);
}
/*-------------------------------------------------------------------------*/
