#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"


/*-------------------------------------------------------------------------*/
void buildA(int DEBUG)
/*-------------------------------------------------------------------------*/
{
    int i;
    double A[length][3];

    for (i=0;i<nobject;i++) 
    {
        if ((object[i])&&(fluidbound[i])) 
        {
	        constructA(A,i);
	        buildLU(A,i);
	        if(DEBUG) debugA(A,i); 
	    }
    }
} 
/*-------------------------------------------------------------------------*/
void constructA(
                double A[length][3], 
                int obj)
/*-------------------------------------------------------------------------*/
{
int j;
double f=0.0;
double delt, invdelt;

invdelt=1.0/delta[0];
delt=delta[0];  /* this is a constant size mesh */

f=0;
if ((MULTI)&&(MIDEAR)) {
	if (incident[obj][SM]) f=density*area[SM][0]/Tm2mm_ms;
	if (incident[obj][ST]) f+=density*area[ST][0]/Tm2mm_ms;
	}
if ((!MULTI)&&(MIDEAR)) {
	if (incident[obj][SM]) f=2.0*density*area[SM][0]/Tm2mm_ms;
	}

A[0][DIAG]=alpha[obj][0]*delt + f + invdelt;  
A[0][SUBDIAG]=(-invdelt)-vgamma[obj][0] * 0.5;

/* this should never be used but its defined anyway, just in case...? */
A[0][SUPRADIAG]=(-invdelt)+vgamma[obj][0] * 0.5;

for (j=1;j<length;j++) { 
	A[j][DIAG]=alpha[obj][j]*delt + 2.0*invdelt;
    A[j][SUPRADIAG]=(-invdelt)-vgamma[obj][j] * 0.5;
	A[j][SUBDIAG]= (-invdelt)+vgamma[obj][j] * 0.5;
    	} 

} 
/*-------------------------------------------------------------------------*/
void buildLU(
             double A[length][3], 
             int obj) 
/* decompose A into its upper (U) and lower (L) matrices                   */
/*-------------------------------------------------------------------------*/
{
int j;

/* construct L matrix */
L[obj][0][LDIAG]=sqrt((double) A[0][DIAG]);
L[obj][1][LSUBDIAG]=A[0][SUPRADIAG]/L[obj][0][LDIAG];

for (j=1;j<(length-1);j++) {
	L[obj][j][LDIAG]=sqrt((double) (A[j][DIAG] - L[obj][j][LSUBDIAG]*L[obj][j][LSUBDIAG]));
    L[obj][j+1][LSUBDIAG]= A[j][SUPRADIAG]/L[obj][j][LDIAG];
    } 
L[obj][length-1][LDIAG]=sqrt((double) (A[length-1][DIAG]- (L[obj][length-1][LSUBDIAG]*L[obj][length-1][LSUBDIAG])));

/* construct U matrix */
for (j=0;j<length;j++) {
    U[obj][j][UDIAG]=L[obj][j][LDIAG];
    U[obj][j][USUPRADIAG]= L[obj][j+1][LSUBDIAG];
    } 

} 
/*-------------------------------------------------------------------------*/
void debugA(
            double A[length][3], 
            int obj)
/*-------------------------------------------------------------------------*/
{
int j;
static FILE *fp;
static open;

if (!open)
	{
	if ((fp = fopen("debugA","w"))==NULL) 
		error("Error opening debug file in debugA()","in A.c");
	open=TRUE;
	}

fprintf(fp,"A matrix\n");
for (j=0;j<length;j++)  fprintf(fp,"|%e\t%e\t%e| vgamma=%e \n",A[j][0],
A[j][1],A[j][2], vgamma[obj][j]); 
fprintf(fp,"\n");

fprintf(fp,"Lower triangular matrix L\n");
for (j=0;j<length;j++)  fprintf(fp,"|%e\t%e| \n",L[obj][j][0],L[obj][j][1]); 
fprintf(fp,"\n");

fprintf(fp,"Upper triangular matrix U\n");
for (j=0;j<length;j++)  fprintf(fp,"|%e\t%e| \n",U[obj][j][0],U[obj][j][1]); 
fprintf(fp,"\n");

}
/*-------------------------------------------------------------------------*/
