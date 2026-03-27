#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"

/*-------------------------------------------------------------------------*/
void netprint(FILE *fp)
/*-------------------------------------------------------------------------*/
{
    int i,j, nobj;
    int stablehi, stablelow;
    int OHCprinted = FALSE;
    double ca1, ca2, cb1, cb2, rb1, rb2;
    double m1, m2, s1, s2, r1, r2, d1, d2, l1, l2; 
    double ss1, ss2, rr1, rr2;
    double hifreq, lofreq;

    fprintf(fp,"\n");
    fprintf(fp,"\n/********************************************************************/");
    fprintf(fp,"\n\tNet object dimensions and parameter ranges");
    fprintf(fp,"\n/********************************************************************/");

    nobj=0;
    for (i=0;i<nobject;i++) if (object[i]) nobj++;

    fprintf(fp,"\n\n\tMesh size:\t\t %d by %d\n",length,nobj);

    fprintf(fp,"\n\tZweig's N: \t\t %e.\n",ZweigN);
  
    for (i=0;i<nobject;i++) 
    {        
	    if ((object[i]) && (i != ELMO)) 
        {
		    m1=M[i][0]; m2=M[i][length-1];
		    s1=S[i][0]; s2=S[i][length-1];
		    r1=R[i][0]; r2=R[i][length-1];
		    d1=damp[i][0]; d2=damp[i][length-1];
		    l1=lambda[i][0]; l2=lambda[i][length-1];

		    hifreq=fpx(s1,m1,r1,&stablehi); 
		    lofreq=fpx(s2,m2,r2,&stablelow);

		    fprintf(fp,"\n\t");
		    nameit(i,fp);
		    fprintf(fp,"\n");
		    fprintf(fp,"\n\t");
	        fprintf(fp, "       base                             apex");
	        fprintf(fp,"\n\t");
	        fprintf(fp, " ----------------------------------------------------");
	        fprintf(fp,"\n\t");
	        fprintf(fp," %e  <=       Mass       <=   %e", m1,m2);
	        fprintf(fp,"\n\t");
	        fprintf(fp," %e  <=     Stiffness    <=   %e", s1,s2);
	        fprintf(fp,"\n\t");
	        fprintf(fp," %e  <=    Resistance    <=   %e", r1,r2);
	        fprintf(fp,"\n\t");
	        fprintf(fp," %e  <=     Damping      <=   %e", d1,d2);
	        fprintf(fp,"\n\t");
	        fprintf(fp," %e  <=      Lambda      <=   %e", l1,l2);

	        fprintf(fp,"\n\t");
            if (stablehi == imaginary) fprintf(fp,"Imaginary");
            else fprintf(fp," %e",hifreq);
	        fprintf(fp,"  <=   Frequency (Hz)  <=  ");
            if (stablelow == imaginary) fprintf(fp,"Imaginary\n");
            else fprintf(fp," %e \n",lofreq);

	        if (fluidbound[i]) 
            {
		        fprintf(fp,"\n\t");
		        fprintf(fp," %e  <=      Alpha	     <=   %e", alpha[i][0],alpha[i][length-1]);
		        fprintf(fp,"\n\t");
		        fprintf(fp," %e  <=  Membrane width  <=   %e", beta[i][0], beta[i][length-1]);
	    	    fprintf(fp,"\n");
            }
            if ((i == OHC) || ((motileOHC)&&(nobj <3)))
            {
                if (!OHCprinted)
                {
                    OHCprinted = TRUE;

                    ca1 = Ca0 * exp(Ca1 * mesh[0]);
                    ca2 = Ca0 * exp(Ca1 * mesh[length-1]);
                    cb1 = Cb0 * exp(Cb1 * mesh[0]);
                    cb2 = Cb0 * exp(Cb1 * mesh[length-1]);
                    rb1 = Rb0 * exp(Rb1 * mesh[0]);
                    rb2 = Rb0 * exp(Rb1 * mesh[length-1]);
	               
                    fprintf(fp,"\n\t");
	                fprintf(fp, " ----------------------------------------------------");

	                fprintf(fp,"\n\t");
	                fprintf(fp," %e  <=        Ca        <=   %e", ca1,ca2);
	                fprintf(fp,"\n\t");
	                fprintf(fp," %e  <=        Cb        <=   %e", cb1,cb2);
	                fprintf(fp,"\n\t");
	                fprintf(fp," %e  <=        Rb        <=   %e", rb1,rb2);
	                fprintf(fp,"\n\t");
	                fprintf(fp," %e  <=    Channel N     <=   %e", stereo_N[0],stereo_N[length-1]);
	                fprintf(fp,"\n\t");
	                fprintf(fp," %e  <=      Length      <=   %e", ohc_l[0],ohc_l[length-1]);
	                fprintf(fp,"\n\t");
	                fprintf(fp," %e  <=   Cilia height   <=   %e", stereo_h[0],stereo_h[length-1]);
	                fprintf(fp,"\n\t");
	                fprintf(fp," %e  <= Pivot stiffness  <=   %e", sp0*exp(sp1*mesh[0]),sp0*exp(sp1*mesh[length-1]));
	                fprintf(fp,"\n\t");
                }
            }
	    }
    }
    if (CUTE) 
    { 
        fprintf(fp,"\n\n");
        for (i=0;i<nobject;i++) 
        { 
            for (j=i;j<nobject;j++) 
            { 
                if ((object[i])&&(coupled[i][j])) 
                {
	                fprintf(fp,"\n\t");
	                nameit(i,fp);
	                fprintf(fp," coupled to ");
	                nameit(j,fp);
	                fprintf(fp,"\n");
	                s1=cS[i][j][0]; s2=cS[i][j][length-1];
	                r1=cR[i][j][0]; r2=cR[i][j][length-1];
	                ss1=cS[j][i][0]; ss2=cS[j][i][length-1];
	                rr1=cR[j][i][0]; rr2=cR[j][i][length-1];
	                if ((ss1!=s1)||(ss2!=s2)||(rr1!=r1)||(rr2!=r2))
		                fprintf(fp, "\n\t 		Coupling is nonsymmetric");
	                printcZ(fp,s1,s2,r1,r2);
                }
            }
        }
    }
    fprintf(fp,"\n\n");
    for (i=0;i<nchannel;i++) 
    {
	    if ((channel[i])&&((i==SM)||(MULTI))) 
        {
		    fprintf(fp,"\n\t");
		    namechannel(i,fp);
		    fprintf(fp,"\n\n\t");
		    fprintf(fp, "       base                             apex");
		    fprintf(fp,"\n\t");
		    fprintf(fp, " ----------------------------------------------------");
		    fprintf(fp,"\n\t");
		    fprintf(fp," %e  <=   Scala area   <=   %e",
		    area[i][0], area[i][length-1]);
		    fprintf(fp,"\n\t");
            if (viscosity) {
		        fprintf(fp," %e  <=     Gamma0      <=   %e",
		        Gamma0[i][0], Gamma0[i][length-1]);
		        fprintf(fp,"\n\t");
		        fprintf(fp," %e  <=     Gamma1      <=   %e",
		        Gamma1[i][0], Gamma1[i][length-1]);
		        fprintf(fp,"\n\t");
		        fprintf(fp," %e  <=     vgamma      <=   %e",
		        vgamma[i][0], vgamma[i][length-1]);
		        fprintf(fp,"\n\n");
            }
	    }
    }
}
/*-------------------------------------------------------------------------*/
void printcZ(
             FILE *fp,
             double s1,
             double s2,
             double r1,
             double r2)
/*-------------------------------------------------------------------------*/
{
	fprintf(fp,"\n\t");
	fprintf(fp, "       base                             apex");
	fprintf(fp,"\n\t");
	fprintf(fp, " ----------------------------------------------------");
    if ((s1 != 0) || (s2 != 0))
    {
	    fprintf(fp,"\n\t");
	    fprintf(fp," %e  <=      Stiffness    <=   %e", s1,s2);
    }
        if ((r1 != 0) || (r2 != 0))
    {
	    fprintf(fp,"\n\t");
	    fprintf(fp," %e  <=     Resistance    <=   %e", r1,r2);
        }
	fprintf(fp,"\n\t");
}
/*-------------------------------------------------------------------------*/
