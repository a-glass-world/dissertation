/**************************************************************************
 *  ohcadmin.c
 **************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <direct.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <io.h>
#include <math.h>

#include "gold.h"

#define ohcdir  "ohcfiles/"
#define ohcglob "ohcglobals"

/*-------------------------------------------------------------------------*/
void tellOHCparms()
/*-------------------------------------------------------------------------*/
{
    struct fileinfo ohcglobal;

    strcpy(ohcglobal.c,ohcglob);
    if ((ohcglobal.fp=fopen(ohcglobal.c,"w"))==NULL) 
    {
	    fprintf(stderr,"file open for %s  failed. \n",ohcglobal.c);
	    error("\tin ohcadmin.c, tellOHCparms()\n","\t ohc file open failed\n");
    }
    ohcglobal.open=TRUE;
    list_ohc(ohcglobal.fp);
    fclose(ohcglobal.fp);
}
/*-------------------------------------------------------------------------*/
void writeOHCfile(
                  double x,
                  double y,
                  char n[],
                  int fn)
/*-------------------------------------------------------------------------*/
{
    
	if (fn >= nohcfile) 		
		error("\tin ohcadmin.c, writeOHCfile()\n","\t indexed out of array.. \n");

	if (!ohcfile[fn].open) 
		openOHCfiles(n,fn);

	fprintf(ohcfile[fn].fp,"%e %e\n",x,y);
}
/*-------------------------------------------------------------------------*/
void closeOHCfiles()
/*-------------------------------------------------------------------------*/
{
    int fn;

    for (fn = 0; fn < nohcfile; fn++) 
    {
	    if (ohcfile[fn].open)
            fclose(ohcfile[fn].fp);
        ohcfile[fn].open = FALSE;
	}
}
/*------------------------------------------------------------------------*/
void openOHCfiles(
                  char n[],
                  int i)
/*-----------------------------------------------------------------------*/
{

	makeohcname(n,i);
	if ((ohcfile[i].fp=fopen(ohcfile[i].c,"w"))==NULL) {
		fprintf(stderr,"file open for %s  failed. \n",ohcfile[i].c);
		error("\tin ohcadmin.c, openOHCfiles()\n",
			"\t ohc file open failed\n");
			}
	ohcfile[i].open=TRUE;
}
/*-----------------------------------------------------------------------*/
void makeohcname(
                  char n[],
                  int i)
/*-----------------------------------------------------------------------*/
{
    char objname[20];

    strcpy(ohcfile[i].c,currentdirname);
    strcpy(objname,"/OHC"); 

    strcat(ohcfile[i].c,objname);
    strcat(ohcfile[i].c,n);
}
/*-----------------------------------------------------------------------*/
