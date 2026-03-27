/**************************************************************************
   trans.c
 **************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"
/* ----------------- used only by elmo and dell and trans ----------------*/
int writej;
extern int writej;

double timelimit;
extern double timelimit;

/*-------------------------------------------------------------------------*/
void capacitivemembrane(
                     double Ra,
                     double x,
                     int YIN,
                     int DOUT,
                     double t,
                     int i)
/*-------------------------------------------------------------------------*/
{
    double f, g;

    /* solve the circuit */
    f = ((Vm * meshRb[i]) + (V0 * Ra)) - (Ra + meshRb[i]) * Y[YIN][0][ELMO].loc[i];
    g =  (Ra * meshRb[i]) * (meshCa0[i] + meshCb0[i]);

    /* set the derivative and scale */
    D[DOUT][0][ELMO].loc[i] = f*1e3/g;
}
/*--------------------------------------------------------------------------*/
void resistivemembrane(
                     double Ra,
                     double x,
                     int YIN,
                     int DOUT,
                     double t,
                     int i)
/*--------------------------------------------------------------------------*/
{

    double f, g, y;


    y = (Ra - restingRa[i])/restingRa[i];

    f = shape_factor_a[i] * y * (Vm - V0);
    g = (1 + shape_factor_a[i]) * (1 + shape_factor_a[i] + y);

    D[DOUT][0][ELMO].loc[i] = f/g;

}
/*--------------------------------------------------------------------------*/
double VbRest(
              double Ra,
              double Rb)
/*--------------------------------------------------------------------------*/
{
    return((Vm*Rb + V0*Ra)/(Ra + Rb));
}
