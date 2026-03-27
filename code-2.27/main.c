/*
 *  main.c
 *
 *  Program to simulate time-domain response of n-dimensional nonlinear 
 *  dynamics of cochlea.  
 *
 *  $ Revision: 1.0  $
 *  $ Date:  30 April 1996 $
 *  $ Author:  M. E. Corey  $
 *  $ Location:  /usrs/albert/gold/main.c $ 
 *	Initial Revision
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <direct.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <io.h>

#include "gold.h"

#define solveit   TRUE /* turn off and on the kernal simulation */
#define slash	  0
#define under	  1


/*------------------------------------------------------------------------*/
void main()
/*------------------------------------------------------------------------*/
{

    initUpperControls();

    mainEventLoop();

    fprintf(stderr,"Bye...\n");
}
/*------------------------------------------------------------------------*/
void initUpperControls()	 		/* let the user choose...*/
/*------------------------------------------------------------------------*/
{
    int quick, isnapt;

    quick=query("Set controls to defaults? No (0). Yes (1).",0,1);

    if (quick) 
    {
	    CUTE        = TRUE;
        ranging     = TRUE;
	    MULTI       = FALSE;
        viscosity   = FALSE;
	    estimating  = FALSE;
	    pauses      = 0;
	    scanning    = FALSE;
	    inputtype   = scanInput;
	    vary        = FALSE;
	    snapt       = 15;
      
	}
    else 
    {
	    CUTE        = query("Uncoupled (0) or coupled? (1)",0,1);
	    MULTI       = query("Regular (0) or multichannel (1) strategies?",0,1);
	    viscosity   = query("Inviscid (0) or viscous (1) fluids?",0,1);
	    ranging     = query("Range through the parameters? No(0), Yes (1)",0,1);
        if (ranging)
        {
            estimating = FALSE;
            scanning = FALSE;
        }
        else 
        {
	        estimating  = query("Constant (0) or estimated (1) coefficients ?",0,1);
	        if (estimating) 
            {
    	        pauses  = query("No pause?  (0) Pause at parm change? (1) Pause at cycle end? (2)",0,2);
	        }
	        else 
            {
                pauses  = 0;
            }

	        scanning    = query("Single run (0) or data scan (1)?",0,1);

        }
	    if (!scanning) 
            inputtype = fixedInput; 
        else 
            inputtype = scanInput;

	    vary    = query("Constant step (0) or variable rk4 (1) ?",0,1); 
	    isnapt  = query("Enter time when snapshot should be taken (will be divided by 10)",1,500);
	    snapt   = ((double) isnapt)/10;
	}
}
/*------------------------------------------------------------------------*/
void mainEventLoop()
/*------------------------------------------------------------------------*/
{
    int done, edone, starting, err;
    int i;
    int newpco = TRUE;
    int dummy;
    /* p is parm, c is coef, o is obj, co is coupled object */
    int p = 0;
    int c = 0;
    int o = 0;
    int rindex =0;

    static int stability;
    static int mtype, ninput, datatype;
    static int matparmdone;
    static int op, oc, oo, obs;  /* keep track of old parameter indices */
    static int obs_is_set;

    /* logical initializations */
    starting    = TRUE;
    edone       = FALSE;
    done        = FALSE;
    stability   = stable;
    err         = FALSE;

    /* index initializations */
    ninput = 0;
    nmodel = 0;

    while (!done) 
    {

	    if (starting) 
        {
		    userinitplan(&datatype,&mtype);
		    systeminitplan(mtype);
		}

	    defineinput(ninput,inputtype,datatype);

        /* open an estimate file for handling overview of scanning results */
        if ((!estimating) && (scanning) && (ninput == 0))
        {
            openestimate(datatype);
        }

	    if (estimating) 
        {
	        if ((ninput==0)||(stability!=stable)) 
            {
	            if (starting)
                {
		            newparm=TRUE;
		            openestimate(datatype);
		            for (i=0;i<nobject;i++) 
		                if (object[i]) shortsummary(estimatefile.fp,i,0);
		        }

	            edone=buildestimate(stability,mtype,&p,&c,&o);

		        if ((p!=op)||(c!=oc)||(o!=oo)||(starting)) 
			        newpco=TRUE; 
		        else 
                    newpco=FALSE;

		        if (!solveit) 
                { /* don't solve it just show the outputs...*/
	    	        summary(stderr,datatype,ninput,nmodel,solveit);
			        shortsummary(stderr,o,p);
	    	    }

		        if (!edone) printMOREstuff(stderr, p, c, o);
		        op=p;oc=c;oo=o;
	        }
        }
        else if (ranging)
        {
            if (starting)
            {
                o=p=c=rindex=0;
            }
            else 
            {
                rindex++;
                if (rindex > 1)
                {
                    rindex = 0;
                    c++;
                }
                if (c>1)
                {
                    c = 0;
                    p++;
                }
                if ((p>2) && (o != nobject) && (!matparmdone))
                {
                    p = 0;
                    o++;
                }
                while ((!matparmdone) && ((!object[o])||(o == ELMO)) && (o < nobject)) o++;
                if (o == nobject)
                {
                    matparmdone = TRUE;
                    o = 0;
                    p = 3;
                    rindex = 0;
                    c = 0;
                }
                if (matparmdone)
                {
                    if (p > 10) 
                        edone = TRUE;
                }  
                /* reset */
                impedance(thesisZ4DOF);
                coupling(thesis4DOFCP);

            }
            if (!edone)
            {
                fprintf(stderr,"--------------------------------------------------------------\n");
                switch(p)
                {
                case MASS:
                    {
                        mass[o][c] = range[rindex][p][o][c];
                        fprintf(stderr,"New mass parameter value is mass[%d][%d] = %e \n",o,c,mass[o][c]);
                        break;
                    }
                case STIFF:
                    {
                        stiff[o][c] = range[rindex][p][o][c];
                        fprintf(stderr,"New stiff parameter value is stiff[%d][%d] = %e \n",o,c,stiff[o][c]);
                        break;
                    }
                case RESIST:
                    {
                        resist[o][c] = range[rindex][p][o][c];
                        fprintf(stderr,"New resist parameter value is resist[%d][%d] = %e \n",o,c,resist[o][c]);
                        break;
                    }
                case cACRLBAMstiff:
                    {
                        cstiff[ACRL][BAM][c] = cstiff[BAM][ACRL][c] = range[rindex][p][nobject][c];
                        fprintf(stderr,"Couple for ACRL and BAM (pillar foot) is cstiff[%d][%d][%d] = %e \n",
                            ACRL,BAM,c,cstiff[ACRL][BAM][c]);
                        break;
                    }
               case cACRLOHCstiff:
                    {
                        cstiff[ACRL][OHC][c] = cstiff[OHC][ACRL][c] = range[rindex][p][nobject][c];
                        fprintf(stderr,"Couple for ACRL and OHC (cuticular) is cstiff[%d][%d][%d] = %e \n",
                            ACRL,OHC,c,cstiff[ACRL][OHC][c]);
                        break;
                    }
               case cACRLTMstiff:
                    {
                        cstiff[ACRL][TM][c] = cstiff[TM][ACRL][c] = range[rindex][p][nobject][c];
                        fprintf(stderr,"Couple for ACRL and TM (cilia) is cstiff[%d][%d][%d] = %e \n",
                            ACRL,TM,c,cstiff[ACRL][TM][c]);
                        break;
                    }
                case cACRLTMresist:
                    {
                        cresist[ACRL][TM][c] = cresist[TM][ACRL][c] = range[rindex][p][nobject][c];
                        fprintf(stderr,"Couple for ACRL and TM (cilia) is cresist[%d][%d][%d] = %e \n",
                            ACRL,TM,c,cresist[ACRL][TM][c]);
                        break;
                    }
               case cACRLACRLstiff:
                    {
                        cstiff[ACRL][ACRL][c] = range[rindex][p][nobject][c];
                        fprintf(stderr,"Couple for ACRL and ACRL (reticular lamina) is cstiff[%d][%d][%d] = %e \n",
                            ACRL,ACRL,c,cstiff[ACRL][ACRL][c]);
                        break;
                    }
               case cOHCOHCstiff:
                    {
                        cstiff[OHC][OHC][c] = range[rindex][p][nobject][c];
                        fprintf(stderr,"Couple for OHC and OHC (phylangeal processes) is cstiff[%d][%d][%d] = %e \n",
                            OHC,OHC,c,cstiff[OHC][OHC][c]);
                        break;
                    }
               case cOHCBAMstiff:
                    {
                        cstiff[OHC][BAM][c] = cstiff[BAM][OHC][c] = range[rindex][p][nobject][c];
                        fprintf(stderr,"Couple for OHC and BAM (Deiters' cell) is cstiff[%d][%d][%d] = %e \n",
                            OHC,BAM,c,cstiff[OHC][BAM][c]);
                        break;
                    }
               case cOHCBAMresist:
                    {
                        cresist[OHC][BAM][c] = cresist[BAM][OHC][c] = range[rindex][p][nobject][c];
                        fprintf(stderr,"Couple for OHC and BAM (Deiters' cell) is cresist[%d][%d][%d] = %e \n",
                            OHC,BAM,c,cresist[OHC][BAM][c]);
                        break;
                    }
                }
            }
        }

	    if (!edone) 
        {
	        /* open the output for the next model */
	        if (estimating)
		        openoutput(ninput,p,c,nmodel,TRUE); 
            else if (ranging)
		        openoutput(o,p,c,rindex,TRUE); 
	        else if (scanning) 
		        openoutput(ninput,nmodel,nmodel,99,TRUE); 
	        else 
		        openoutput(ninput,nmodel,nmodel,99,FALSE); 


	        /* set up an array index for the first (calibration) runs */
            /*  -- obs is an index into arrays of length fftlength -- */
	        if ((estimating) && (calibrating))
                obs=(int) (((double) (3.0/cl))*fftlength);
            else if ((!estimating) && (!obs_is_set))
            {
                obs = setobs(freq[0],mtype);
                obs_is_set = TRUE;
            }

            peakindex = (length/fftlength)*obs;

            /* describe the basic parameters and strategies */
	        planprint(explain.fp,ninput,mtype);

	        if (solveit) 
            { 
		        /* run the simulation */
	            stability = model(starting,obs,mtype,nmodel,c); 

		        /* show the analysis */
		        if (stability==stable)  
                { 
	                    printFund(nmodel);
	                    printFFT(nmodel); 
		        }
            
		        /* build the response structs and success measures */
	            if ((estimating) || (scanning))
                {
		            buildresponse(ninput,p,c,o,newpco,stability,obs);
		            buildsuccess(ninput,datatype,stability);
		            /* I don't even want to see the bad maps */
		            if ((stability == stable) || (stability == unstable)) 
                    {
			            summary(stderr,datatype,ninput,nmodel,solveit);
		    	        summary(estimatefile.fp,datatype,ninput,nmodel,solveit);
                        if (estimating) 
                        {
		    	            shortsummary(stderr,o,p);
		    	            shortsummary(estimatefile.fp,o,p);
                        }
			        }
		        }

		        /* NOW, if you have just finished the calibration run,
		           you have the observation point */
	            if ((estimating) && (calibrating)) 
		            obs = response[ninput][nmodel][BAM].maxindex; 

		        /* and you can say where it is...*/
	            if (stability == stable) comment(ninput);
            }
            else
            {
                /* just set it up and talk about it 
                stability = buildnet(mtype,starting);
                netprint(explain.fp);
                initTD(mtype);
                strategyprint(explain.fp);
                shortsummary(stderr,o,0);
	            shortsummary(stderr,o,3);
	            shortsummary(stderr,o,6); */

            }

	    }

	    /* are we really done? */
	    done=decide(mtype, &ninput, edone,stability,p,c,o);

	    /* yes, we've definitely started... */
	    starting=FALSE;
	}

    if ((estimating) || (scanning)) 
    { 
        /* synopsis and close up */
	    shortsummary(stderr,o,0);
	    shortsummary(stderr,o,3);
	    shortsummary(stderr,o,6);
	    shortsummary(estimatefile.fp,o,0);
	    shortsummary(estimatefile.fp,o,3);
	    shortsummary(estimatefile.fp,o,6);
	    printpfile(mtype);
	    fclose(estimatefile.fp); 
	    if (!ranging) 
            fclose(rrr.fp); /* yep, the brrr, rrr global file */
	}
    else
    {
        /* just close everything */
        closeoldmodel();
    }

}
/*------------------------------------------------------------------------*/
void systeminitplan(int modelnumber)
/*------------------------------------------------------------------------*/
{
    int geo,imp,me,cp;

	/* now do the rest */
	activity(modelnumber);
	topology(modelnumber,&geo,&imp,&me,&cp);
	geometry(geo);
	impedance(imp);
	middleear(me);
	coupling(cp);

}
/*------------------------------------------------------------------------*/
void userinitplan(
                  int *datatype,
                  int *modelnumber)
/*------------------------------------------------------------------------*/
{
	int errcode,whatever;
	unsigned short mode;

	estimateactive = constantactive = FALSE;

	if (estimating) 
		estimateactive=query("Estimate passive (0) or active (1)?",0,1);
	else
	    constantactive=query("Run passive (0) or active (1)?",0,1);

	fprintf(stderr,"What model do you want to run?\n");
	if ((!estimateactive)&&(!constantactive)) 
    {
	    fprintf(stderr,"\t viergever:     passive 1DOF \t\t 0\n");
	    fprintf(stderr,"\t passive1DOF:   passive 1DOF \t\t 1\n");
	    fprintf(stderr,"\t passive2DOF:   passive 2DOF \t\t 2\n");
        fprintf(stderr,"\t passive4DOF:   passive 4DOF \t\t 3\n");
        fprintf(stderr,"\t neely2DOF:     passive 2DOF \t\t 21\n");
	    *modelnumber=query("Please input model number",0,21);
	}
	else 
    {
	    fprintf(stderr,"\t active2DOF: active 2DOF \t\t\t 4\n");
	    fprintf(stderr,"\t active3DOF: active 3DOF \t\t\t 5\n\n");
        fprintf(stderr,"\t motile OHC, capacitive circuit                6\n");
	    fprintf(stderr,"\t motile OHC, extracelluar circuit              7\n");
	    fprintf(stderr,"\t dynamic coupled OHC, capacitive circuit       8\n");
	    fprintf(stderr,"\t dynamic coupled OHC, extracelluar circuit     9\n");
	    fprintf(stderr,"\t active stiffness OHC soma                    10\n");
	    fprintf(stderr,"\t active stiffness OHC cilia                   11\n");
	    fprintf(stderr,"\t active stiffness OHC, both soma and cilia    12\n");
	    fprintf(stderr,"\t motile OHC, somatic stiffness                13\n");
	    fprintf(stderr,"\t motile OHC, cilia stiffness                  14\n");
	    fprintf(stderr,"\t motile OHC, soma and cilia                   15\n");
	    fprintf(stderr,"\t dynamic couple, motile OHC, soma and cilia   16\n");
	    *modelnumber=query("Please input model number",4,16);
	}

	if (estimating) 
    {
	    strcpy(outputdir,datadir); 
	    switch(*modelnumber) 
        {
		    case viergever1DOF: {
			    strcat(outputdir, "/Viergever"); 
			    break; }
		    case thesis1DOF: {
			    strcat(outputdir, "/m1DOF"); 
			    break; }
		    case passive2DOF: {
			    strcat(outputdir, "/pass2DOF"); 
			    break;}
		    case neely2DOF: {
			    strcat(outputdir, "/neely2DOF"); 
			    break;}
		    case active2DOF: {
			    strcat(outputdir, "/act2DOF"); 
			    break;}
		    case active3DOF: {
			    strcat(outputdir, "/act3DOF"); 
			    break;}
		    case passive4DOF: {
			    strcat(outputdir, "/pass4DOF"); 
			    break;}
		    case motile4DOFcapacitive: { 
			    strcat(outputdir, "/motile4DOFcap"); 
			    break;}
		    case motile4DOFextracellular: { 
			    strcat(outputdir, "/motile4DOFextra"); 
			    break;}
		    case couple4DOFcapacitive: { 
			    strcat(outputdir, "/couple4DOFcap"); 
			    break;}
		    case couple4DOFextracellular: { 
			    strcat(outputdir, "/couple4DOFextra"); 
			    break;}
		    case stiff4DOFsoma: { 
			    strcat(outputdir, "/stiff4DOFsoma"); 
			    break;}
		    case stiff4DOFcilia: { 
			    strcat(outputdir, "/stiff4DOFcilia"); 
			    break;}
		    case stiff4DOFboth: { 
			    strcat(outputdir, "/stiff4DOFboth"); 
			    break;}
		    case motilestiff4DOFsoma: { 
			    strcat(outputdir, "/motilestiff4DOFsoma"); 
			    break;}
		    case motilestiff4DOFcilia: { 
			    strcat(outputdir, "/motilestiff4DOFcilia"); 
			    break;}
		    case motilestiff4DOFboth: { 
			    strcat(outputdir, "/motilestiff4DOFboth"); 
			    break;}
            case all4DOF: { 
			    strcat(outputdir, "/all4DOF"); 
			    break;}
		    default: { error("case undefined","initplan()"); }
		}
	} 
	else if (scanning)
    {
		strcpy(outputdir,datadir);
		strcat(outputdir, "/scan"); 
	}
	else if (ranging)
    {
		strcpy(outputdir,datadir);
		strcat(outputdir, "/perturb"); 
	}
	else
    {
		strcpy(outputdir,datadir);
		strcat(outputdir, "/misc"); 
	}

	errcode=_mkdir(outputdir);
	if ((errcode == (-1))&&(errno==EACCES)) 
    {
		fprintf(stderr,"The directory %s already exists.\n",outputdir);
		whatever=query("Quit? No (0), Yes (1)",0,1);
		if (whatever==1) 
            exit(0);
	}
	else if ((errcode == (-1))&&(errno==ENOENT)) 
    {
		fprintf(stderr,"The path to the directory %s was not found.\n",outputdir);
        error("initplan() can not access",outputdir);
	}

	/* set mode to 777 octal */
	mode=7*64 + 7*8 + 7;
	errcode=_chmod(outputdir,_S_IWRITE);

	if (errcode==(-1)) 
    {
            error("main.c, initplan(), chmod can not access",outputdir);
    }

	if (scanning) {
	    if ((estimateactive)||(constantactive)) 
        {
	        fprintf(stderr,"What animal state do you want?\n"); 
	        fprintf(stderr,"\t Live       \t\t 1\n");
	        fprintf(stderr,"\t Furosemide \t\t 2\n");
	        fprintf(stderr,"\t Exposure   \t\t 3\n");
	        animal_status=query("Please input state number",1,4)-1;
		}
	    else 
        {
	        fprintf(stderr,"What animal state do you want?\n"); 
	        fprintf(stderr,"\t Post-mortem\t\t 1\n");
	        fprintf(stderr,"\t Furosemide \t\t 2\n");
	        if (query("Please input state number",1,2) == 1)
			animal_status = post;
	        else
			animal_status = furo;
		}

	    fprintf(stderr,"What type of input do you wish to match?\n");

	    switch(animal_status)
        { 
		    case live: 
            {
			    fprintf(stderr,"\t Click data \t\t\t 0\n");
			    fprintf(stderr,"\t Isointensity data \t\t 1\n");
			    fprintf(stderr,"\t Two-tone distortion data \t 2\n");
			    fprintf(stderr,"\t Two-tone suppression data \t 3\n");
			    fprintf(stderr,"\t Tuning curves \t\t\t 4\n");
			    fprintf(stderr,"\t Quick tune \t\t\t 5\n");
			    fprintf(stderr,"\t Gain curves \t\t\t 6\n");
			    *datatype=query("Please input data type",0,6);
			    break;
			}
		    case loud: 
            {
			    fprintf(stderr,"\t Click data \t\t 7\n");
			    fprintf(stderr,"\t Isointensity data \t 8\n");
			    *datatype=query("Please input data type",7,8);
			    break;
			}
		    case post: 
            {
			    fprintf(stderr,"\t Click data \t\t 9\n");
			    fprintf(stderr,"\t Tuning curves \t\t 10\n");
			    *datatype=query("Please input data type",9,10);
			    break;
			}
		    case furo: 
            {
			    fprintf(stderr,"\t Click data \t\t 11\n");
			    fprintf(stderr,"\t I/O curves \t\t 12\n");
			    fprintf(stderr,"\t Click FFT  \t\t 13\n");
			    *datatype=query("Please input data type",11,13);
			    break;
			}
		    default:
            { 
                error("status case undefined","initplan()");
		    }
	    }

	    if ((*datatype==6)||(*datatype==8)||(*datatype==0)||(*datatype==10)) 
	    {
	         RECIO=query("Input pure click (0) or Reccio's click (1)?",0,1);
	    }

	    appendDname(outputdir,*datatype,slash,currentdir);

    }	
    appendDname(outputdir,DONTCARE,DONTCARE,currentdir);

}
/*------------------------------------------------------------------------*/
void appendDname(
                 char putdir[nchar],
                 int dtype,
                 int addpunc,
                 char getdir[nchar])
/*------------------------------------------------------------------------*/
{

    strcpy(getdir,putdir);

    switch(addpunc) {
	    case slash: 
        {
		    strcat(getdir,"/");
		    break; 
        }
	    case under: 
        {
		    strcat(getdir,"_");
		    break; 
        } 
	    case DONTCARE: 
        {
            break; 
        }
	}
		
    switch(dtype) {
	case 0: case 7: case 9: case 11: { strcat(getdir,"click"); break; }
	case 1: case 8: { strcat(getdir,"isoint"); break; }
	case 2: { strcat(getdir,"dist2tone"); break; }
	case 3: { strcat(getdir,"supp2tone"); break; }
    case 4: case 5: case 10: {strcat(getdir,"tune"); break; }
	case 6: { strcat(getdir,"gain"); break; }
	case 12: { strcat(getdir,"clickFFT"); break; }
	case 13: { strcat(getdir,"io_curve"); break; }
    case DONTCARE: {break;}
	default:{ error("In main.c, in appendDname()","case undefined");}
	}
} 
/*------------------------------------------------------------------------*/
int model(
          int starting,
          int obs,
          int modeltype,
          int nm,
          int coef)
/*------------------------------------------------------------------------*/
{
    int stability;
    static int startedsolve;

    /* discretize the model */
    stability = buildnet(modeltype,starting);
    netprint(explain.fp);

    /* make sure it makes sense, kinda */
    if (estimating) 
        constraints(&stability,obs);
   
    if (stability!=stable) 
        return(stability);
    if (!startedsolve) 
    {
        initTD();
        startedsolve=TRUE;
    }
    else 
    {
        resetTD();
    }

    strategyprint(explain.fp);


    /* SOLVE * SOLVE * SOLVE * SOLVE */
    stability=solve(TRUE,nm); 

    if (stability==unstable) 
    {
        fprintf(stderr,"This model is unstable.\n");
        if (estimatefile.open)
            fprintf(estimatefile.fp,"This model is unstable.\n");
        if (explain.open)
            fprintf(explain.fp,"This model is unstable.\n");
        }

    return(stability);
}
/*------------------------------------------------------------------------*/
int decide(
           int modeltype, 
           int *ninput, 
           int edone,
           int stability,
           int p,
           int c,
           int o)
/*------------------------------------------------------------------------*/
{
    int done;

    done=FALSE;
    if (estimating) 
    {
        if (edone) 
            done=TRUE;
        if (nestimates<=nmodel) 
        {
	        resetparm(p,c,o);
	        done=TRUE;
	        edone=TRUE;
	        fprintf(stderr,"Oops, ran out of memory for estimates");
	    }
        else 
        {
	        if (scanning) 
            {
   	            if ((*ninput<nscanmodel-1)&&(stability==stable)) 
                {
                    (*ninput)++; 
                }
   	            else 
                {
	                (*ninput)=0;
	                (nmodel)++;
	            }
	        }
	        else 
            {
                (nmodel)++;
            }
	    }
    }
    else if (scanning) 
    {
   	    if (*ninput < nscanmodel-1) 
            (*ninput)++; 
	    else 
            done=TRUE;
	}
    else if (ranging)
    {
        if (edone)
            done = TRUE;
    }
    else 
    {
        done=TRUE;
    }

    return(done);
}
/*------------------------------------------------------------------------*/
void openestimate(int dtype)
/*------------------------------------------------------------------------*/
{

    strcpy(estimatefile.c,outputdir);

    if (estimating)
        strcat(estimatefile.c,"/estimate");
    else if (scanning)
        strcat(estimatefile.c,"/scanning");

    appendDname(estimatefile.c,dtype,under,estimatefile.c);

    if ((estimatefile.fp=fopen(estimatefile.c,"w"))==NULL) 
    {
	    fprintf(stderr,"file open for %s  failed. \n",estimatefile.c);
	    error("\tin main()\n","\t estimate file open failed\n");
	}

    estimatefile.open=TRUE;

}
/*------------------------------------------------------------------------*/
void sortresponses(
                   int totmod,
                   int datatype)
/*------------------------------------------------------------------------*/
{
    int i,j,k,x;
    double fit[nestimates];
    double hold;
    struct scanstruct holds[maxscan];
    struct scanstruct holdr[maxscan][nobject];

    for (i=0;i<nscanmodel;i++) 
    { /* totmod is number of models */
        for (j=0;j<totmod;j++) 
        { /* nscanmodel is number of inputs */
	        for (k=0;k<response[i][j][BAM].nfeature;k++) 
            { 
	            fit[i]+=response[i][j][BAM].diff[k]; 
            }
        }
    }

    for (i=0;i<=totmod;i++) 
    {
        for (j=i;j<=totmod;j++) 
        {
	        if (fit[j]<fit[i]) 
            {
	            hold=fit[i];
	            fit[i]=fit[j];
	            fit[j]=hold;
	            for (x=0;x<nscanmodel;x++) 
                {
	    	        holds[x]=success[x][i];
	    	        success[x][i]=success[x][j];
	    	        success[x][j]=holds[x]; 
		            for (k=0;k<nobject;k++) 
                    { 
                        if (object[k]) 
                        {
	    	                holdr[x][k]=response[x][i][k];
	    	                response[x][i][k]=response[x][j][k];
	    	                response[x][j][k]=holdr[x][k]; 
		                } 
                    }
                }
		    }
        }
    }
	     
    for (i=0;i<=5;i++) 
    { 
        for (j=0;j<nscanmodel;j++) 
        { 
	        summary(estimatefile.fp,datatype,j,i,solveit);  
	        summary(stderr,datatype,j,i,solveit); 
	    } 
    }
}
/*------------------------------------------------------------------------*/
void constraints(
                 int *stability,
                 int obs)
/*------------------------------------------------------------------------*/
{

    if (*stability==imaginary) 
    {
        fprintf(stderr,"This model has an imaginary frequency map.\n");
	    fprintf(explain.fp,"This model has an imaginary frequency map.\n");
   	    if (estimating) fprintf(estimatefile.fp,"This model has an imaginary frequency map.\n");
	    *stability = stable; /* its ok */
    }
    else 
    {
	    if ((animal_status==live)&&(estimateactive)) 
        {
		    if ((lomap>1.0)||(himap<15)||(himap>40)) 
			    *stability=badmap;
	    }
    }

    if (*stability==badmap) 
    {
	    fprintf(stderr,"Frequency map is %f < %f\n",lomap, himap);
	    fprintf(estimatefile.fp,"Frequency map is %f < %f\n",lomap,himap); 
    	    fprintf(stderr,"This model has an unusable frequency map.\n");
   	    fprintf(estimatefile.fp,"This model has an unusable frequency map.\n");
    	    fprintf(explain.fp,"This model has an unusable frequency map.\n");
	}

    /* check the observation place for grave errors */
    if ((!calibrating)&&(freqdiff>1.0)&&(estimating)) 
    {
	    fprintf(stderr,"Frequency map is %f < %f\n",lomap, himap);
	    fprintf(stderr,"Bad observation place: %f mm\n",obs*cl/fftlength);
	    fprintf(estimatefile.fp,"Bad observation place: %f mm\n",obs*cl/fftlength);
	    fprintf(explain.fp,"Bad observation place: %f mm\n",obs*cl/fftlength);
	    *stability=badmap;
	}
}
/*------------------------------------------------------------------------*/
void printpfile(int modeltype)
/*------------------------------------------------------------------------*/
{
    struct fileinfo pfile;
    FILE *fp;

    strcpy(pfile.c,outputdir);
    strcat(pfile.c,"/parm"); 

    if ((pfile.fp=fopen(pfile.c,"w"))==NULL) 
    {
	    fprintf(stderr,"file open for %s  failed. \n",pfile.c);
	    error("\tin main()\n","\t print parameter file open failed\n"); 
    }
    else 
    {
	    fprintf(stderr,"file open for %s  succeeded. \n",pfile.c);
	}

    pfile.open=TRUE;
    fp=pfile.fp;

    if (object[BAM])     printparm(fp,BAM);
    if (object[TM])      printparm(fp,TM);
    if (object[OHC]) 
    {
        printparm(fp,OHC);
        printOHCparm(fp);
    }

    if (object[ACRL])    printparm(fp,ACRL);

    if (coupled[BAM][TM])    printcparm(fp,TM,BAM);
    if (coupled[BAM][OHC])   printcparm(fp,BAM,OHC);
    if (coupled[BAM][ACRL])  printcparm(fp,ACRL,BAM);
	if (coupled[TM][ACRL])   printcparm(fp,ACRL,TM);
	if (coupled[OHC][ACRL])  printcparm(fp,OHC,ACRL);
	if (coupled[OHC][OHC])   printcparm(fp,OHC,OHC);
	if (coupled[ACRL][ACRL])   printcparm(fp,ACRL,ACRL);

    fclose(pfile.fp);
}
/*------------------------------------------------------------------------*/
void printparm(FILE *fp,
               int k)
/*------------------------------------------------------------------------*/
{
	fprintf(fp,"%e\n",mass[k][0]);
	fprintf(fp,"%e\n",mass[k][1]);
	fprintf(fp,"%e\n",stiff[k][0]);
	fprintf(fp,"%e\n",stiff[k][1]);
	fprintf(fp,"%e\n",resist[k][0]);
	fprintf(fp,"%e\n",resist[k][1]);
}
/*------------------------------------------------------------------------*/
void printcparm(
                FILE *fp,
                int obj,
                int cobj)
/*------------------------------------------------------------------------*/
{
	fprintf(fp,"%e\n",cmass[obj][cobj][0]);
	fprintf(fp,"%e\n",cmass[obj][cobj][1]);
	fprintf(fp,"%e\n",cstiff[obj][cobj][0]);
	fprintf(fp,"%e\n",cstiff[obj][cobj][1]);
	fprintf(fp,"%e\n",cresist[obj][cobj][0]);
	fprintf(fp,"%e\n",cresist[obj][cobj][1]);
}
/*------------------------------------------------------------------------*/
void printOHCparm(FILE *fp)
/*------------------------------------------------------------------------*/
{
	fprintf(fp,"%e\n",sp0);
	fprintf(fp,"%e\n",sa);
	fprintf(fp,"%e\n",gc);
	fprintf(fp,"%e\n",sl);
	fprintf(fp,"%e\n",stereo_N0);
	fprintf(fp,"%e\n",stereo_N1);
}
/*------------------------------------------------------------------------*/
int setobs(double infreq, int mtype)
/*------------------------------------------------------------------------*/
{
/* for reference:
    passive system, 4DOF, 500 Hz =~ 30 */

    infreq *= 1000; /* convert to Hz (hence, an integer) */
    switch((int) infreq)        
    {
        case 100: /* 100 Hz */
        {
            switch(mtype) 
            {
                case thesis1DOF: return(100);
                default: {}
            } 
            break;
        }

        case 1000: /* 1.000 KHz */
        {
            switch(mtype) 
            {
                case thesis1DOF: return(56);
                default: {}
            } 
            break;
        }

        case 10000: /* 10.000 KHz */
        {
            switch(mtype) 
            {
                case thesis1DOF: return(12);
                default: {}
            } 
            break;
        }

    }
   
    fprintf(stderr,"For frequency %f, fftlength %d, the observation index is undefined\n",
        infreq,fftlength);
    return(query("Please input the observation index",0,fftlength));

}
/*------------------------------------------------------------------------*/
