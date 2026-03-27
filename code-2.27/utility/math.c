
/*------------------------------------------------------------------------*/
/*  				math.c 					*/
/*------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"

/*------------------------------------------------------------------------*/
double exp10(double power)
/*------------------------------------------------------------------------*/
{
double ten, value;

ten=10.0;
value=pow(ten,power);
return(value);
}
/*------------------------------------------------------------------------*/
struct corey_complex cmult(
                           struct corey_complex a, 
                           struct corey_complex b)
/*------------------------------------------------------------------------*/
{
struct corey_complex c;

c.re=(a.re*b.re)-(a.im*b.im);
c.im=(a.im*b.re)+(a.re*b.im);

return(c);
}
/*------------------------------------------------------------------------*/
struct corey_complex cdiv(
                          struct corey_complex a, 
                          struct corey_complex b)
/*------------------------------------------------------------------------*/
{
struct corey_complex c, cnvB, denom, num;

cnvB.re=b.re;
cnvB.im=(-b.im); 

denom=cmult(b,cnvB);
num=cmult(a,cnvB);

c.re=num.re/denom.re;
c.im=num.im/denom.re;


return(c);
}
