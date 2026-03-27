/**************************************************************************
 *  summary.c
 **************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"

#define datap		1
#define responsep 	2	
#define successp	3	
#define allp

/*------------------------------------------------------------------------*/
void quickfile(
               int dtype,
               int inp,
               int mod)
/*------------------------------------------------------------------------*/
{
    static int started;
    static int lastmodel;
    int j;
    struct scanstruct x, y, z;

    x = data[inp];
    y = response[inp][mod][BAM];
    z = success[inp][mod];

    /* make sure the response file is open */
    if (!started) 
    {
        switch(dtype) 
        {
	    case postclick:    { newfile(mod,"BAM_","postclick",FALSE);  break; }
	    case posttune:     { newfile(mod,"BAM_","posttune",FALSE);   break;  }
	    case livetune:     { newfile(mod,"BAM_","livetune",FALSE);   break;  }
	    case quicklivetune:{ newfile(mod,"BAM_","quicktune",FALSE);  break;  }
	    case liveclick:    { newfile(mod,"BAM_","liveclick",FALSE);  break; }
	    case isointensity: { newfile(mod,"BAM_","iso",FALSE);        break; } 
	    default: 
            {
	        fprintf(stderr,"Warning:  generic global response files \n");
	        newfile(mod, "BAM", "", FALSE);
	        }
	    }
        started = TRUE;
    } 

    /* print the data and response */
    if (!y.cal) 
    { 
        if (mod != lastmodel) 
        {
    	    fprintf(rrr.fp,"Estimate %d, ",mod);    
            switch(y.estparm) 
            {
	        case RESIST: 
                { 
	            fprintf(rrr.fp,"\t resist[%d][%d] = %e (estimated %e) \n",
	            y.estobj,y.estcoef,resist[y.estobj][y.estcoef],y.est); 
	            break; 
                }
	        case STIFF: 
                { 
	            fprintf(rrr.fp,"\t stiff[%d][%d]  = %e (estimated %e) \n",
	            y.estobj,y.estcoef,stiff[y.estobj][y.estcoef],y.est); 
	            break; 
                }
	        case MASS: 
                { 
	            fprintf(rrr.fp,"\t mass[%d][%d]   = %e (estimated %e) \n",
	            y.estobj,y.estcoef,mass[y.estobj][y.estcoef],y.est); 
	            break; 
                }
	        default: {}
	        }
        }
        for (j=0;j<y.nfeature;j++) 
        {
	        fprintf(rrr.fp,"\t\t| %f ",y.feature[j]); 
	        fprintf(rrr.fp," - %f = ",x.feature[j]); 
	        fprintf(rrr.fp," %f |\n",z.feature[j]); 
	    }
        fprintf(rrr.fp,"\n");
    } 

}
/*------------------------------------------------------------------------*/
void newfile(
             int num,
             char text1[],
             char text2[],
             int showmodelnumber)
/*------------------------------------------------------------------------*/
{
    char n[4],blank[2];

    /* close it if its open */
    if (rrr.open) fclose(rrr.fp);

    /* clobber it if there's anything in there */
    blank[0]='\0';
    strcpy(rrr.c,blank);

    /* put the file directory name in */
    strcpy(rrr.c,globaldir);

    /* add the text given by the parameter list */
    blank[0]='/';
    blank[1]='\0';
    strcat(rrr.c,blank);
    strcat(rrr.c,text1);
    strcat(rrr.c,text2);

    /* add the number given by the parameter list */
    n[0]=num/100 + '0';
    n[1]=num/10-num/100 + '0';
    n[2]=num%10 + '0';
    n[3]='\0';
    if (showmodelnumber) strcat(rrr.c,n); 

    /* tell what the file name is and open it */
    if ((rrr.fp = fopen(rrr.c,"w")) == NULL) 
    { 
        fprintf(stderr,"file name %s \n",rrr.c);
        error("in newfile() \n","Open failed for global data file.\n");
    }
    /* set logical open to true */
    rrr.open=TRUE;

}
/*------------------------------------------------------------------------*/
void summary(
             FILE *fp,
             int dtype,
             int inp,
             int mod,
             int solved)
/*------------------------------------------------------------------------*/
{
    int i,k,p,c;
    int co1, co2;
    double obsloc, maxloc,norm;
    struct scanstruct x,y,z;

    x=data[inp];
    y=response[inp][mod][BAM];
    z=success[inp][mod];

    if (solved) 
    {
        fprintf(fp,"\n--------------------------------------------------------------------\n");
	    fprintf(fp,"Estimate %d, Datatype %d, Input %d is %e dB and %f KHz\n",
	    mod,dtype,x.inputnum,x.dB,x.freq);
	    if (x.cal) fprintf(fp,"\t This is a calibration run.\n");
    }
    else 
    {
	    if (inp==0) fprintf(fp,"Estimate %d, Datatype %d\n", mod,dtype);
    }

    if (solved) 
    {
        obsloc=(((double) y.obsindex)/fftlength)*cl;
        fprintf(fp,"\t Observation location is \t\t %e (index %d)\n",obsloc,y.obsindex);

        maxloc=(((double) y.maxindex)/fftlength)*cl;
        fprintf(fp,"\t CP                      \t\t %e\n",maxloc);
        if (estimating)
        {
            switch(y.estparm) 
            {
	            case RESIST:    { fprintf(fp,"\t Damping"); break; }
	            case STIFF:     { fprintf(fp,"\t Modular stiffness"); break; }
	            case MASS:      { fprintf(fp,"\t Plate thickness"); break;} 
	            default:        { fprintf(fp,"\t %s",ui_returnname(y.estparm)); break;}
	        }
            if (y.estparm > 2) 
            {   
                getcoupledobjs(y.estparm,&co1,&co2);
                fprintf(fp,"[%d][%d]         \t\t %e \n", y.estobj, y.estcoef, y.est); 
            }
            else
            {
                fprintf(fp,"[%d][%d][%d]         \t\t %e \n", co1,co2, y.estcoef, y.est); 
            }
        }

        for (k = 0; k < x.nfeature; k++) 
        {
	        fprintf(fp,"\t feature %d \n",k);
            fprintf(fp,"\t\t data                    \t %e \n",x.feature[k]);
            fprintf(fp,"\t\t simulated response      \t %e \n",y.feature[k]);
	        fprintf(fp,"\t\t data - response         \t %e \n",z.feature[k]);
	        if (estimating)
            {
                fprintf(fp,"\t\t delta success (ds)      \t %e \n",z.diff[k]);
                fprintf(fp,"\t\t delta estimate (de)     \t %e \n",z.intval);
                fprintf(fp,"\t\t derivative (ds/de)      \t %e \n",z.deriv[k]);
            }
	    }

        if (!x.cal) 
        {
            /*
	        fprintf(fp,"\t Success (normed)\t\t\t %e\n",z.norm);
	        */
	        if (x.inputnum == nscanmodel-1) 
            {
	            fprintf(fp,"\n");
	            p = y.estparm;
	            c = y.estcoef;

	            /* skip the first one, it doesn't count */
    	        for (k = 1; k < nscanmodel; k++) 
                {
		            if (sensitive.nmod[p][c] != 0)
	                    fprintf(fp,"\t Sensitive (input %d) \t\t\t %e\n",k,
                                sensitive.featinp[p][c][k]/sensitive.nmod[p][c]);
		        }

    	        for (k = 0; k < z.nfeature; k++) 
                {
                    if (sensitive.nmod[p][c] != 0)
	                    fprintf(fp,"\t Sensitive (feature %d)\t\t\t %e\n",k,
                                sensitive.feature[p][c][k]);
                }

	            norm=0.0;
    	        for (k = 1; k < nscanmodel; k++) 
	    	        for (i = 0; i < data[inp].nfeature; i++) 
		                norm += success[k][mod].feature[i]/data[inp].nfeature;
	            if ((nscanmodel - 1) > 0) norm /= (nscanmodel-1);

	            fprintf(fp,"\t Frequency-place map\t\t\t %f to %f\n",lomap,himap);
	            fprintf(fp,"\t Success (Features) \t\t\t %f\n",norm);
	            fprintf(fp,"\t Success (Obs is %f)\t\t\t %f\n",obsloc,freqdiff);
	            fprintf(fp,"\t Success (Total) \t\t\t %e\n",z.scan_norm_tot);
	            fprintf(fp,"\n");
	        }
	    }
        /* now write the summary files */
        if (fp == stderr) 
            quickfile(dtype,inp,mod);
    } 
}
/*------------------------------------------------------------------------*/
void summaryscan(
                 FILE *fp,
                 struct scanstruct x,
                 int mod,
                 int dtype,
                 int flag)
/*------------------------------------------------------------------------*/
{
    int k;
    double obsloc, maxloc;

    fprintf(fp,"Estimate %d, Datatype %d, Input %d is %e dB and %f KHz\n",
	    mod,dtype,x.inputnum,x.dB,x.freq);
    if (x.cal) 
	    fprintf(fp,"\t This is a calibration run.\n");

    obsloc=(((double) x.obsindex)/fftlength)*cl;
    fprintf(fp,"\t Observation location is \t\t %e (index %d)\n",obsloc,x.obsindex);
    fprintf(fp,"\t CF of observation place \t\t %e\n",x.observe);
    if ((flag==responsep)||(flag==successp)) 
    {
        maxloc=(((double) x.maxindex)/fftlength)*cl;
        fprintf(fp,"\t CP                      \t\t %e\n",maxloc);

        if (x.newpm) 
            fprintf(fp,"\n\t New parameter type.\n");
        else 
            fprintf(fp,"\n\t Continuing with same parameter type.\n");

        switch(x.estparm) 
        {
	        case RESIST:{fprintf(fp,"\t Resist[%d][%d]         \t\t\t %e \n",
	                      x.estobj,x.estcoef,x.est); break; }
	        case STIFF: {fprintf(fp,"\t Stiff[%d][%d]          \t\t\t %e \n",
	                      x.estobj,x.estcoef,x.est); break; }
	        case MASS:  {fprintf(fp,"\t Mass[%d][%d]           \t\t\t %e \n",
		             x.estobj,x.estcoef,x.est); break; }
	        default: {}
	    }
        fprintf(fp,"\t Estimate interval       \t\t %e\n",x.intval);
    }
    for (k=0;k<x.nfeature;k++) 
    {
        if (flag!=successp) 
        { 
            fprintf(fp,"\t feature[%d]             \t\t %e \n",k, x.feature[k]);
	    }
        else 
        {
	    fprintf(fp,"\t diff=response - data[%d]\t\t\t %e \n",k,x.feature[k]);
            fprintf(fp,"\t diff[%d]                \t\t %e \n",k,x.deriv[k]);
            fprintf(fp,"\t deriv[%d]               \t\t %e \n",k,x.deriv[k]);
	    }
    }
    if (flag==successp) 
    {
        fprintf(fp,"\t normed differences -feat \t\t %e\n",x.norm);
        fprintf(fp,"\t normed differences -feat & input\t %e\n",x.scan_norm_tot);
    } 
}
/*------------------------------------------------------------------------*/
void shortsummary(
                  FILE *fp,
                  int obj,
                  int parmtype)
/*------------------------------------------------------------------------*/
{
    int i,j;
    int co1, co2;

    if (parmtype < 3) 
    {
        nameit(obj,fp);
        fprintf(fp," material parameters\n"); 
        for (j=0;j<ncoef;j++) 
        { 
	        fprintf(fp,"\t"); 
	        fprintf(fp,"resist[%d][%d]=%e",obj,j,resist[obj][j]); 
	        fprintf(fp,"\n");
	    }
        for (j=0;j<ncoef;j++) 
        { 
	        fprintf(fp,"\t"); 
	        fprintf(fp,"stiffness[%d][%d]=%e",obj,j,stiff[obj][j]); 
	        fprintf(fp,"\n");
	    }
        for (j=0;j<ncoef;j++) 
        { 
	        fprintf(fp,"\t"); 
	        fprintf(fp,"mass[%d][%d]=%e",obj,j,mass[obj][j]); 
	        fprintf(fp,"\n");
	    }
    }

    if (parmtype > 2) 
    {
        getcoupledobjs(parmtype,&co1,&co2);
        fprintf(fp,"Couple for "); nameit(co1,fp); fprintf(fp," and "); nameit(co2,fp);
        fprintf(fp,"\n");
        for (j=0;j<ncoef;j++) 
	        fprintf(fp,"\t cresist[%d][%d][%d] = %e\n ",co1,co2,j,cresist[co1][co2][j]); 
        for (j=0;j<ncoef;j++) 
	        fprintf(fp,"\t cstiff[%d][%d][%d] = %e\n ",co1,co2,j,cstiff[co1][co2][j]); 
    }
        
    fprintf(fp,"\n");
    return;
}
/*------------------------------------------------------------------------*/
