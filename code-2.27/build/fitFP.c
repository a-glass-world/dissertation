#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"

/*------------------------------------------------------------------------*/
void fitFP(
           double low,
           double high,
           double apex,
           int obj)
/*------------------------------------------------------------------------*/
{
double r1,r2,m1,m2,k1,k2;
double wlow, whigh, wbase2, wapex2, m,r,k; 
double l, checklow,checkhigh;
double twopi;
int dummy;

twopi=M_PI*2.0;
/* adjust to radians */
wlow=low*twopi;
whigh=high*twopi;
wbase2=whigh*whigh;
wapex2=wlow*wlow;

r1=resist[obj][0];
r2=resist[obj][1];

m1=mass[obj][0];
m2=mass[obj][1];

if ((m1==0)&&(m2==0)) 
	error("in fitFP.c, subroutine fitFP\n","Mass undefined \n");

k1=(wbase2*4.0*m1*m1 + r1*r1)/(4.0*m1);
m=m1*exp(m2*apex);
r=r1*exp(r2*apex);
l=(wapex2*4.0*m*m +r*r)/(4.0*m*k1);
k2=log(l)/apex;
k=k1*exp(k2*apex);

dummy=0;
checkhigh=fpx(k1,m1,r1,&dummy);
checklow=fpx(k,m,r,&dummy);

stiff[obj][0]=k1;
stiff[obj][1]=k2;
}
/*------------------------------------------------------------------------*/
double fpx(
           double s,
           double m,
           double r,
           int *stabi)
/*------------------------------------------------------------------------*/
{
    double result;
    double dv;
    double a,b,c,d,e;
    double twopi;

    /* (4ms-r^2)/4m^2 */

    twopi=M_PI*2.0;

    dv  = twopi;
    a   = 4.0 * m * s; 
    b   = r * r;
    c   = 4.0 * m * m;
    d   = a - b;
    e   = d/c;

    if (e<0.0) 
    {
	    result = 0;
	    *stabi = imaginary;
    }
    else 
    {
        result = 1000*sqrt(e)/dv;
        *stabi = stable;
    }

    return(result);
}
/*--------------------------------------------------------------------------*/

