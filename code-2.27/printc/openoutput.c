#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <direct.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <io.h>
#include <math.h>

#include "gold.h"


/*------------------------------------------------------------------------*/
void openoutput(
                int n,
                int m,
                int c,
                int nm,
                int aut)
/*------------------------------------------------------------------------*/
{
    char modeldirname[nchar];
    char modelnum[5];

    closeoldmodel();

    makemodeldir(n,m,c,nm,modeldirname,aut);

    if (!estimating) 
	    fprintf(stderr,"The model directory name is %s \n",modeldirname);

    strcat(modeldirname,explainname);

    if (nm != DONTCARE)
    {
        modelnum[0]='_';
        modelnum[1]=nm/100+'0';
        modelnum[2]=(nm%100)/10+'0';
        modelnum[3]=nm%10+'0';
        modelnum[4]='\0';
    }
    else
    {
        modelnum[0]='\0';
    }

    strcat(modeldirname,modelnum);
    strcpy(explain.c,modeldirname);

    if ((explain.fp=fopen(explain.c,"w"))==NULL) 
    {
        fprintf(stderr,"Open failed for file %s\n",explain.c);
	    error("in openoutput.c \n","openoutput().\n");
    }

    explain.open=TRUE;

}
/*------------------------------------------------------------------------*/
void makemodeldir(
                  int n,
                  int nm, 
                  int c,
                  int rindex,
                  char modeldirname[nchar],
                  int aut)
/*------------------------------------------------------------------------*/
{
    char m[5], inname[100];
    char blank[1];
    int errcode, whatever;
    unsigned short mode;
    static int started;

    mode=7*64 + 7*8 +7;

    blank[0]='\0';
    strcpy(modeldirname,blank);

    if (aut) 
    {
        strcpy(modeldirname,currentdir);
        if (!ranging)
        {
            m[0]='/';
            m[1]=n/10+'0';
            m[2]=n%10+'0';
            m[3]='\0';
            strcat(modeldirname,m);
        }
        else
        {
            if (rindex == 0)
                strcat(modeldirname,"/lo_");
            else if (rindex == 1)
                strcat(modeldirname,"/hi_");
            else
                error("more garbage","openoutput.c");

            if (nm < cACRLBAMstiff)
            {
                switch(n)
                {
                    case BAM:{strcat(modeldirname,"BAM");break; }
                    case OHC:{strcat(modeldirname,"OHC");break;}
                    case ACRL:{strcat(modeldirname,"ACRL");break;}
                    case TM:{strcat(modeldirname,"TM");break;}
                    case nobject:{break;}
                    default:
                    {
                        error("no, no, no, what object is this??","openoutput.c");
                    }
                }
            }
        }
        switch(nm) 
        {
	        case 0: {strcat(modeldirname,"_R");;break;}
	        case 1: {strcat(modeldirname,"_S");;break;}
	        case 2: {strcat(modeldirname,"_M");;break;}
	        case 3: {strcat(modeldirname,"_ACRL_BAM_S");break;}
	        case 4: {strcat(modeldirname,"_ACRL_OHC_S");break;}
	        case 5: {strcat(modeldirname,"_ACRL_TM_R");break;}
	        case 6: {strcat(modeldirname,"_ACRL_TM_S");break;}
	        case 7: {strcat(modeldirname,"_ACRL_ACRL_S");break;}
	        case 8: {strcat(modeldirname,"_OHC_OHC_S");break;}
	        case 9: {strcat(modeldirname,"_OHC_BAM_R");break;}
	        case 10: {strcat(modeldirname,"_OHC_BAM_S");break;}
	        default:{error("Oops, finish case", "openoutput.c"); break;}
	    }
        switch(c) 
        {
	        case 0: {strcat(modeldirname,"0"); break;}
	        case 1: {strcat(modeldirname,"1"); break;}
	        default:{}
	    }
    }
    else
    {
        fprintf(stderr,"Please input a directory run for output \n>> ");
        scanf("%s",inname);
        if ((inname[0]=='q') && (inname[1]=='\0')) exit(0);
        strcat(modeldirname,currentdir);
        m[0]='/';
        m[1]='\0';
        strcat(modeldirname,m);
        strcat(modeldirname,inname);
    }

    errcode=_mkdir(modeldirname);
    if ((errcode!=0)&&(!started)) 
    { 
	    fprintf(stderr,"%s already exists.\n",modeldirname);
	    fprintf(stderr,"Continuing may overlay files in this directory.\n");
	    whatever=query("Quit? No (0), Yes (1)",0,1);
	    if (whatever==1) 
            exit(0);
	}
    started=TRUE;

    /* set mode to 777 octal */
    errcode=_chmod(modeldirname,_S_IWRITE);
    if (errcode!=0) 
	    error("ERROR: openoutput.c, makemodeldir()\n\tcan not access",
        modeldirname);
    /* otherwise, save the directory name for posterity */
    strcpy(currentdirname,modeldirname);
}
/*------------------------------------------------------------------------*/
void closeoldmodel()
/*------------------------------------------------------------------------*/
{
int i,j;
if (explain.open) {
	fclose(explain.fp);
	explain.open=FALSE;	
	}

for (i=0;i<nfile;i++) {
	for (j=0;j<actmem;j++) {
		if (fftfile[i][j].open) {
			fclose(fftfile[i][j].fp);
			fftfile[i][j].open=FALSE;	
			}
		if (fundfile[i][j].open) {
			fclose(fundfile[i][j].fp);
			fundfile[i][j].open=FALSE;	
			}
		}  
	}
for (i=0;i<nobject;i++) {
	for (j=0;j<ndebug;j++) {
		if (debug[i][j].open) {
			fclose(debug[i][j].fp);
			debug[i][j].open=FALSE;	
	} } }

for (i=0;i<nobject;i++) {
	if (snapx[i].open) {
		fclose(snapx[i].fp);
		snapx[i].open=FALSE;	
		}
	if (snapv[i].open) {
		fclose(snapv[i].fp);
		snapv[i].open=FALSE;	
		}
	if (clickx[i].open) {
		fclose(clickx[i].fp);
		clickx[i].open=FALSE;	
		} 
	if (clickv[i].open) {
		fclose(clickv[i].fp);
		clickv[i].open=FALSE;	
		} 
	if (clickt[i].open) {
		fclose(clickt[i].fp);
		clickt[i].open=FALSE;	
		} 
	}
  closeOHCfiles();
}
/*------------------------------------------------------------------------*/
