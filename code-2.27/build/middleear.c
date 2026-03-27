#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"
/*------------------------------------------------------------------------------*/
void middleear(int MEtype)
/*------------------------------------------------------------------------------*/
{

	switch(MEtype) {
		case neelyME: {
			MIDEAR=TRUE;
			Tm2mm_ms= 23.55;
			Tm= 2.0;
			T2m= 4.0;
			Am= 34.7;
			Rm=14.6;
			Sm= 1817.0;
			T2mRm=T2m*Rm;
			break;
			}
		case diepenME: {
			MIDEAR=TRUE;
			Tm2mm_ms= 7.92e-01;
			Tm= 2.0;
			T2m= 4.0;
			Am= 60.0;	/* mm*mm */
			Rm=6.22;	/* resistance in mg/ms */
			Sm= 125.08;	/* stiffness in mg/(ms*ms) */
			T2mRm=T2m*Rm;
			break;
			}
		case noME: { 
			MIDEAR=FALSE;
			break;
			}
		default:{
			error("in middleear.c, middleear()","Middle ear transfer function undefined");
			}
		}
}
/*------------------------------------------------------------------------------*/

