/*-------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"

#define matext	".mat"
/* -------------------------------------------------------------------------*/
void fluiddebug(
                double time,
                int x, 
                int y,
                int OB)
/* -------------------------------------------------------------------------*/
{
int i;
char ntype[4];
int inK;
static int count;

if (count > 200) 
{
    error("\tThat's all, folks!\n","\tExiting prematurely\n");
}

count++;


if (OB!=BAM) return;

for (i=0; i< ndebug; i++) {
   if (!debug[OB][i].open) {
	strcpy(debug[OB][i].c,currentdirname);
	strcat(debug[OB][i].c,debugname);
	ntype[0]=OB+'0';	
	ntype[1]=i+'0';
	ntype[2]='\0';	
	strcat(debug[OB][i].c,ntype);
	strcat(debug[OB][i].c,matext);
	if ((debug[OB][i].fp=fopen(debug[OB][i].c,"w"))==NULL)
		error("in deBugPrint.c \n","Open failed for debug file.\n");
	else 
		fprintf(stderr,"The debug files are in %s \n",debug[OB][i].c);
	debug[OB][i].open=TRUE;
	}
    }

fprintf(debug[OB][0].fp,"%e\n",time) ;
fprintf(debug[OB][1].fp,"%e\n",totsignal) ;
fprintf(debug[OB][2].fp,"%e\n",Gme) ;

inK=10;
fprintf(debug[OB][3].fp," %e \n ",G[inK]);
fprintf(debug[OB][4].fp," %e \n ",Gc[inK]);
fprintf(debug[OB][5].fp," %e \n ",PHI[inK]);
fprintf(debug[OB][6].fp," %d \n ",inK);
fprintf(debug[OB][7].fp," %e \n ",time);

}
    
/*-------------------------------------------------------------------------*/
