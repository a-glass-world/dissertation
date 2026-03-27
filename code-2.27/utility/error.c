/*************************************************************************
 *  Errors.c
 **************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"

/*--------------------------------------------------------------------------*/
void error(
           char message1[200],
           char message2[200])
/*--------------------------------------------------------------------------*/
{

if (explain.open) fprintf(explain.fp,"%s %s",message1,message2);
fprintf(stderr,"\n\tMeep, meep.... \n Exiting with an error ...\n\t");
fprintf(stderr," \n %s %s \n",message1,message2);
if (estimating) {
	if (estimatefile.open) fclose(estimatefile.fp);
	if (rrr.open) fclose(rrr.fp);
	}
/* close everything else */
closeoldmodel();

exit(0); 

}	
/*--------------------------------------------------------------------------*/

