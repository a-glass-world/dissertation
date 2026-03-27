#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"

void error();

/*----------------------------------------------------------------------------*/
void snapprint(
               int YIN,
               double time,
               int inloc,
               int st,
               int nm, 
               int newx)
/*----------------------------------------------------------------------------*/
{
    int i,j;


    if ((st == SNAP) && (!estimating) &&(!scanning) && (!ranging)) 
    {
        for (i = 0;i < nobject;i++) 
        { 
            if ((object[i]) && (i != ELMO))
            {
                if (newx) 
                { 
	            snapx[i]=opensnapfiles(nm,i,"x_"); 
	            snapv[i]=opensnapfiles(nm,i,"v_"); 
	            }
	            /* this is the response across the length of the cochlea
                   NOW (at one time) */
	            for (j=0;j<length;j++) {
		            fprintf(snapx[i].fp,"%e %e \n", (double) j*cl/length, Y[YIN][0][i].loc[j]);
		            fprintf(snapv[i].fp,"%e %e \n", (double) j*cl/length, Y[YIN][1][i].loc[j]);
	            }
	        }
        }
    }
    else if (st == CLICK)
    { /* this is the response at one location over a period of time */
        for (i=0;i<nobject;i++) 
        { 
            if ((object[i]) && (i != ELMO))
            {
                if (newx) 
                {
	                if ((!estimating) && (!ranging))
                        clickx[i]=opensnapfiles(nm,i,"x."); 
                    clickv[i]=opensnapfiles(nm,i,"v."); 
	            }
	            if ((!estimating) && (!ranging))
                    fprintf(clickx[i].fp,"%e %e \n", time, Y[YIN][0][i].loc[inloc]);
	            fprintf(clickv[i].fp,"%e %e \n", time, Y[YIN][1][i].loc[inloc]);
            }
        }
    }
}
/*----------------------------------------------------------------------------*/
struct fileinfo opensnapfiles(
                              int nm,
                              int obj,
                              char it[10])
/*----------------------------------------------------------------------------*/
{
    char objname[10];
    char mod[4];
    struct fileinfo ifi;

    ifi.open=0;
    ifi.c[0]='\0';

    if (estimating) 
    {
        mod[0]=nm/100+'0';
        mod[1]=nm/10-nm/100+'0';
        mod[2]=nm%10+'0';
        mod[3]='\0';
    }
    else
    {
        mod[0]='\0';
    }

    switch(obj) 
    {
	    case BAM: { strcpy(objname,"BAM"); break; }
	    case TM: {  strcpy(objname,"TM"); break; }
	    case OHC: { strcpy(objname,"OHCcentroid"); break; }
	    case ACRL: { strcpy(objname,"ACRL"); break; }
	    case ELMO: { strcpy(objname,"ELMO"); break; }
	    default:{error("case undefined in printc/snapprint.c",
		    "snapfiles not open"); }
    }
    strcpy(ifi.c,currentdirname);
    strcat(ifi.c,"/");
    strcat(ifi.c,it);
    strcat(ifi.c,objname);
    strcat(ifi.c,mod);
    if ((ifi.fp=fopen(ifi.c,"w"))==NULL) {
	    fprintf(stderr, "File name %s\n", ifi.c);
	    error("in snapprint.c \n","Open failed.\n");
	    }
    ifi.open=TRUE;

    return(ifi);
}
/*----------------------------------------------------------------------------*/
