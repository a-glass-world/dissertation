#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"

/* before you try to unify fftprint and componentprint, notice they are
   printing different files when estimating/scanning and when not */
#define magnitude   0

#define db          1
#define phaze       2
#define freq_comp   3 
#define imagine     4
#define reall       5


/*-----------------------------------------------------------------------*/
void printFFT(int m)
/*-----------------------------------------------------------------------*/
{
    int n,i,j,k;
    int cont;
    int fundindex;
    int numfile;

    if (!tds.print3Dfft)  return;

    if ((estimating) || (scanning) || (ranging))
    {
        numfile = 1;
    }
    else
    {
        numfile = 4;
    }

    for (i = 0; i < numfile; i++) 
    {
	    for (j = 0; j < actmem; j++) 
        {
		    if (!fftfile[i][j].open) 
            {
			    cont = makename(i, j, m);

			    if (cont) 
                {
				    if ((fftfile[i][j].fp = fopen(fftfile[i][j].c,"w")) == NULL)
                    {
					    error("\tin printfft.c,printFFT()\n","\t fft file open failed\n");
				    }
                    fftfile[i][j].open = TRUE;
			    }
		    }
	    }
    }

    for (n = 0; n < nobject; n++)  
    {
	    if (objfft[n] < actmem) 
        { 
		    i = objfft[n];  
            getfourierindex(&fundindex, n, 1.0);

		    for (j = 0; j < fftlength; j += 2) 
            {
			    for (k = 0; k < TIMESTEPS/2; k++)  
                {
                    gnuprintvalue(i, j, k, fundindex); 
                }
			    for (k = 0; k < nfile; k++) 
                {
                    if (fftfile[k][i].open) 
                    {
                        fprintf(fftfile[k][i].fp,"\n");
				    }
			    }
		    }
	    }
    }
}

/*-----------------------------------------------------------------------*/
void gnuprintvalue(
                int i,
                int j,
                int k,
                int fundindex)
/*-----------------------------------------------------------------------*/
{
    double zero = 0.0;
    double loc = (double) j*cl/fftlength;
    double ffreq;
    if (fundindex != 0)
        ffreq = (((double) k)/((double) fundindex)) * freq[0];
    else
        ffreq = k;

    fprintf(fftfile[magnitude][i].fp,"%e %e %e \n", loc, ffreq, fftFreq[i][j].mag[k]);

    if ((!estimating) && (!scanning) && (!ranging))
    {
        fprintf(fftfile[phaze][i].fp,"%e %e %e \n", loc, ffreq, fftFreq[i][j].phase[k]);

        if (fftFreq[i][j].mag[k]>0.0) 
	        fprintf(fftfile[db][i].fp,"%e %e %e \n", loc, ffreq, 20.0*log10(fftFreq[i][j].mag[k]));
        else 
            fprintf(fftfile[db][i].fp,"%e %e %e \n", loc, ffreq, zero);

        fprintf(fftfile[freq_comp][i].fp,"%e %d %e \n",  loc, k, ffreq);
    }

}
/*-----------------------------------------------------------------------*/
int makename(
             int i,
             int j,
             int m)
/*-----------------------------------------------------------------------*/
{
char num[3];
char objname[20];
char modnum[5];
int k, cont = FALSE;


	for (k = 0; k < nobject; k++) 
    {
        if (objfft[k] == j) 
        {
            getobjectname(objname,k);
            cont = TRUE;
		}
    }
	if (!cont) return(FALSE);

    if (estimating)
    {
	    modnum[0] = m/100 + '0';
	    modnum[1] = m/10-m/100 + '0';
	    modnum[2] = m%10 + '0';
	    modnum[3] = '\0';
    }
    else 
    {
        modnum[0] = '\0';
    }

    strcpy(fftfile[i][j].c,currentdirname);
    strcat(fftfile[i][j].c,fftname);
    strcat(fftfile[i][j].c,objname);
    strcat(fftfile[i][j].c,modnum);

    getextension(num,i);
	
	strcat(fftfile[i][j].c,num);
	return(TRUE);
}
/*----------------------------------------------------------------------------*/
