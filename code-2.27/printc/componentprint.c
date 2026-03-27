#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"

#define magnitude   0
#define db          1
#define phaze       2
#define freq_comp   3 

#define imagine     4
#define reall       5

/*-----------------------------------------------------------------------*/
void printFund(
               int m)
/*-----------------------------------------------------------------------*/
{
    int n,i,j,k;
    int numfile, cont;

    if ((!estimating) && (!scanning) && (!ranging))
        numfile = nfile; 
    else 
        numfile = 3;


    for (i = 0; i < numfile; i++) 
    {
        for (j = 0; j < actmem; j++) 
        {
	        if (!fundfile[i][j].open) 
            {
	            cont = makefundname(i,j,m);
	            if (cont) 
                {
	                if ((fundfile[i][j].fp=fopen(fundfile[i][j].c,"w"))==NULL)
		            {
		                fprintf(stderr,"File name %s \n",fundfile[i][j].c);
		                error("\tin componentprint.c,printFund()\n",
			              "\t fundamental file open failed\n");
			        }
		            fundfile[i][j].open=TRUE;
		        }
	        }
	    }
    }

    for (n = 0; n < nobject; n++) 
    { 
        if (objfft[n] < actmem) 
        { 
	        i = objfft[n];  

	        for (j = 0; j < fftlength; j++) 
            {
	            for (k = 0;k <TIMESTEPS/2.0 ;k++) 
                { 
	                if (fftFreq[i][j].freq[k] == 1.0) 
                    {
		                printfundvalue(i, j, k); 
		            }
		        } 
	            for (k = 0; k < numfile; k++) 
                {
                    fprintf(fundfile[k][i].fp,"\n");
	            }
	        }
        }
    }
}
/*-----------------------------------------------------------------------*/
void printfundvalue(
                    int i,
                    int j,
                    int k)
/*-----------------------------------------------------------------------*/
{
    double zero = 0;

    /* db */
    if (fftFreq[i][j].mag[k] > 0.0) 
	    fprintf(fundfile[db][i].fp,"%e %e ",(double) j*cl/fftlength, 20.0*log10(fftFreq[i][j].mag[k]));
    else 
        fprintf(fundfile[db][i].fp,"%e %e ",(double) j*cl/fftlength, zero);

    /* phase */
    fprintf(fundfile[phaze][i].fp,"%e %e ", (double) j*cl/fftlength, fftFreq[i][j].phase[k]);
    fprintf(fundfile[magnitude][i].fp,"%e %e ", (double) j*cl/fftlength, fftFreq[i][j].mag[k]);

    if ((!estimating) && (!scanning) && (!ranging))
    {
	    fprintf(fundfile[reall][i].fp,"%e %e ", (double) j*cl/fftlength, fftFreq[i][j].real[k]);
	    fprintf(fundfile[imagine][i].fp,"%e %e ", (double) j*cl/fftlength, fftFreq[i][j].im[k]);
	    fprintf(fundfile[freq_comp][i].fp,"%e %e ", (double) j*cl/fftlength, fftFreq[i][j].freq[k]);
	}


}
/*-----------------------------------------------------------------------*/
int makefundname(
                 int i,
                 int j,
                 int m)
/*-----------------------------------------------------------------------*/
{
char num[5];
char objname[20];
char modnum[5];
int k, cont = FALSE;

	cont = FALSE;
	for (k = 0;k < nobject; k++)
    {
        if (objfft[k] == j) 
        {
		    cont=TRUE;
            getobjectname(objname,k);
		}
    }

	if (!cont) return(FALSE); 

    if (estimating) 
    {
	    modnum[0] = '_';
	    modnum[1] = m/100 + '0';
	    modnum[2] = m/10  -  m/100 + '0';
	    modnum[3] = m%10  + '0';
       	modnum[4] = '\0';
    }
    else
    {
       	modnum[0]='\0';
    }

    strcpy(fundfile[i][j].c, currentdirname);
    strcat(fundfile[i][j].c, fundname);
    strcat(fundfile[i][j].c, objname);
    strcat(fundfile[i][j].c, modnum);

    getextension(num,i);
	
	strcat(fundfile[i][j].c,num);
	return(TRUE);
}
/*----------------------------------------------------------------------------*/
void getobjectname(
             char objname[20],
             int obj)

/*----------------------------------------------------------------------------*/
{
    switch(obj) 
    {
	    case BAM: 
            {
                strcpy(objname,"BAM"); 
                break;
            }
		case TM:  
            {
                strcpy(objname,"TM"); 
                break;
            }
		case OHC: 
            {
                strcpy(objname,"OHC"); 
                break;
            }
		case ACRL: 
            {
                strcpy(objname,"ACRL"); 
                break;
            }
		case ELMO: 
            {
                strcpy(objname,"ELMO"); 
                break;
            }
		case OHCSTEREO: 
            {
                strcpy(objname,"OHCSTEREO"); 
                break;
            }
        case DELL: 
            {
                strcpy(objname,"DELL"); 
                break;
            }
		default:
            {
                error("\tin printfft.c, makename()\n","\tundefined objectname case\n");
			    break; 
            }
    }
}
/*----------------------------------------------------------------------------*/
void getextension(
             char num[5],
             int exttype)

/*----------------------------------------------------------------------------*/
{
    
	num[0] = '.';
	switch(exttype) 
    {
		case db: 
            {num[1]='d'; break;}
		case phaze: 
            {num[1]='p'; break;}
		case magnitude: 
            {num[1]='m'; break;}
        case imagine: 
            {num[1]='i'; break;}
		case reall: 
            {num[1]='r'; break;}
        case freq_comp: 
            {num[1]='f'; break;}
		default:
        {
            fprintf(stderr,"Extension type %d",exttype);
            error("\tin componentprint.c, getextension()\n",
                "\tundefined case\n");
            break; 
        }
	}

	num[2]='\0';
}
