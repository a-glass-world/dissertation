/**************************************************************************
 *  buildsuccess.c
 **************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"

#define debugging FALSE

#define Some_cycle_number_never_reached 9999



/*------------------------------------------------------------------------*/
void buildsuccess(
                  int ninput,
                  int dtype,
                  int stability)
/*------------------------------------------------------------------------*/
{
int err;
int nweight;
struct topology intersect, partitions;

	err = makesuccessfeature(ninput,nmodel);
    if (err) return;
	setsensitive(ninput,nmodel);
	measuresuccess(ninput,nmodel,dtype,nweight);

    /* stubs */
	if (ninput == (nscanmodel - 1)) 
        measuretrajectory(nmodel);
	hyperplane(intersect);
	hessian(intersect,partitions);
	bifurcate(partitions);
}

/*------------------------------------------------------------------------*/
int makesuccessfeature(
                        int n,
                        int m)
/*------------------------------------------------------------------------*/
{
    double x,y,in;
    int i;
    int err = FALSE;

    /* you need information from response */
    success[n][m] = response[n][m][BAM];

    /* if the features aren't calculated because this is a calibration or an unstable model, get out */
    if ((response[n][m][BAM].cal)||(response[n][m][BAM].stability != stable))
    {
        if (response[n][m][BAM].cal) 
        {
	        /* all peaks should be about 3.0 mm from base */
	        x = (double) success[0][m].maxindex;
	        y = (double) ((x * cl)/ (double) fftlength);
	        freqdiff = sqrt((y - 3.0) * (y - 3.0));
        }
        err = TRUE;

        success[n][m].nfeature      = 0;
    	success[n][m].scan_norm_tot = 0;
    	success[n][m].norm          = 0; 

        for (i = 0;i < data[n].nfeature; i++) 
        { 
            success[n][m].feature[i] = 0.0;
        }
        return(err);
    }


    for (i = 0;i < data[n].nfeature; i++) 
    { 
	    x = data[n].feature[i];
	    y = response[n][m][BAM].feature[i];
        /* find the difference */
	    in = fabs(x - y);
        /* norm the difference */
	    success[n][m].feature[i] = fabs(in/x);
    }
    return(err);
}
/*------------------------------------------------------------------------*/
void resetsensitivity(
                      int p,
                      int c,
                      int n,
                      int m)
/*------------------------------------------------------------------------*/
{
int i;
	/* reset everything and return */
    	for (i = 0; i < data[n].nfeature; i++) 
    	{
        	success[n][m].deriv[i] = 0.0;	
        	success[n][m].diff[i]  = 0.0;
	}

    /* this says sensitivity has never been calculated */
	if (sense_cycle[p][c] == Some_cycle_number_never_reached)
	{
	    for (i = 0; i<nscanmodel; i++) 
		    sensitive.featinp[p][c][i] = 0.0; 

        for (i = 0; i < data[n].nfeature; i++) 
		    sensitive.feature[p][c][i] = 0.0; 
	}

}
/*------------------------------------------------------------------------*/
int calculable(
               int n,
               int m)
               /* n is input number, m is model number */
/*------------------------------------------------------------------------*/
{

    /* this is not an estimated model */
    if (!estimating)
        return(FALSE);
    /* last simulation was a calibration */
    if (response[n][m][BAM].cal) 
	    return(FALSE);
    /* last simulation was a calibration */
    if (n == 0) 
	    return(FALSE);
    /* this is the first simulation */
    if (m == 0) 
	    return(FALSE);
    /* this is a new parameter */
    if (response[n][m-1][BAM].estcoef != response[n][m][BAM].estcoef)
	    return(FALSE);
    /* this is a new parameter */
    if (response[n][m-1][BAM].estparm != response[n][m][BAM].estparm)
	    return(FALSE);
    /* the last model wasn't stable */
    if (response[0][m-1][BAM].stability!=stable)
	    return(FALSE);

    return(TRUE);
}
/*------------------------------------------------------------------------*/
int sensitivity_calculated(
                           int p,
                           int c)
/*------------------------------------------------------------------------*/
{
    return(cycle > sense_cycle[p][c]);
}
/*------------------------------------------------------------------------*/
void setsensitive(
                  int n,
                  int m)
/*------------------------------------------------------------------------*/
{
    int i,j,k;
    int p,c;

	/* find out what is being estimated */
	p = response[n][m][BAM].estparm;
	c = response[n][m][BAM].estcoef;

	/* for the first model, set it up so we know that the sensitivity has not been calculated */
	if (m == 0)
	    for (i = 0; i < nparm; i++)
		    for (j = 0; j < ncoef; j++)
		        sense_cycle[i][j] = Some_cycle_number_never_reached;


	/* check if we've already done this, if so, get out */
	if (sensitivity_calculated(p,c))
		return;

    /* if we can't do this, reset things */
	if (!calculable(n,m)) 
	{
		resetsensitivity(p,c,n,m);
    	return;
	}


	/* calculate the derivatives */
	for (k=0;k<data[n].nfeature;k++) 
	{
    	success[n][m].diff[k]  = fabs(success[n][m].feature[k] - success[n][m-1].feature[k]);
        if (success[n][m].intval != 0.0)
    	    success[n][m].deriv[k] = success[n][m].diff[k]/success[n][m].intval;
        else
            success[n][m].deriv[k] = 0.0;

	}

    /* count the number of models we measure sensitivity for - 
      this gives a divisor for featinp and tells us whether or not to print */
	if (n == 1) 
		sensitive.nmod[p][c]++;

	/* calculate the sensitivities for initial mesh values */
    for (i = 0; i < data[n].nfeature; i++) 
    {
        /* add all feature differences to input sensitivity */ 
    	sensitive.featinp[p][c][n] += success[n][m].diff[i];
    	/* add each feature difference to its feature sensitivity */ 
    	sensitive.feature[p][c][i] += success[n][m].diff[i];
    }

	sense_cycle[p][c] = cycle;
}
/*------------------------------------------------------------------------*/
void measuresuccess(
                    int n,
                    int m,
                    int dtype,
                    int nweight)
/*------------------------------------------------------------------------*/
{
    int i;

    for (i = 0; i < data[n].nfeature; i++) 
    { 
	    success[n][m].norm += success[n][m].feature[i]/data[n].nfeature; 
	}

    if ((n == nscanmodel-1) && (n > 0))
    {
        for (i = 1; i < nscanmodel; i++) 
        { 
	        success[n][m].scan_norm_tot += success[i][m].norm/(nscanmodel - 1);
	    }
    }
}

/*------------------------------------------------------------------------*/
void measuretrajectory(int m)
/*------------------------------------------------------------------------*/
{
}
/*------------------------------------------------------------------------*/
void hyperplane(struct topology intersect)
/*------------------------------------------------------------------------*/
{
}
/*------------------------------------------------------------------------*/
void hessian(
            struct topology intersect,
            struct topology partitions)
/*------------------------------------------------------------------------*/
{
}
/*------------------------------------------------------------------------*/
void bifurcate(struct topology partitions)
/*------------------------------------------------------------------------*/
{
}
/*------------------------------------------------------------------------*/
