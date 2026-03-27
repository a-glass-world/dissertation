/**************************************************************************
 *  readclick.c
 **************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gold.h>

/* these are strictly local to this file */
#define filelength 1000

struct fileinfo clickfile;
extern struct fileinfo clickfile;

double timex[filelength], magnit[filelength];
extern double timex[filelength], magnit[filelength];
/*------------------------------------------------------------------------*/
double getclick(double ntime)
/*------------------------------------------------------------------------*/
{
int i;
static int linecount;

    if (clickfile.open!=TRUE) {
	 openclick();
	 linecount=readclick(clickfile.fp);
	 }

    for (i=1;i<linecount;i++) {	
        if ((timex[i-1]<=ntime)&&(ntime<timex[i])) {
		return(magnit[i-1]);
		}
	}
return(0.0);
}
/*------------------------------------------------------------------------*/
void openclick()
/*------------------------------------------------------------------------*/
{
strcpy(clickfile.c,clickdir);
strcat(clickfile.c,"/click.txt");

if ((clickfile.fp=fopen(clickfile.c,"r"))==NULL) {
	fprintf(stderr,"file open for %s  failed. \n",clickfile.c);
	error("\tin main()\n","\t click file open failed\n");
	}
clickfile.open=TRUE;
}
/*------------------------------------------------------------------------*/
int readclick(FILE *fp)
/*------------------------------------------------------------------------*/
{
int rval, count;
double x,y;

rval=0;
count=0;
while (rval!=(-1)) {
    rval = fscanf(clickfile.fp,"%lf %lf",&x,&y);
    timex[count] = x;
    magnit[count]= y;
    count++;
    }

return(count-1);
}
/*------------------------------------------------------------------------*/
