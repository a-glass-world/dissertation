#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"

/*----------------------------------------------------------------------------*/
void topology(
              int modeltype,
              int *geo,
              int *imp,
              int *me,
              int *cp)
/*----------------------------------------------------------------------------*/
{
    int i,j, memorylimit;

    memorylimit = 99;
    actmem = 1;

    for (i=0;i<nobject;i++) 
    {
	    object[i]=FALSE;
	    objfft[i]=memorylimit;
	    fluidbound[i]=FALSE;
	    for (j=0;j<nchannel;j++)  incident[i][j]=FALSE;
	    for (j=0;j<nobject;j++) coupled[i][j]=FALSE;
    }

    for (i=0;i<nchannel;i++)  
        channel[i]=FALSE;

    /* sign conventions */
	    UC[BAM]=SM;
	    LC[BAM]=ST;
	    UC[TM]=SM;
	    LC[TM]=SS;


    switch(modeltype) 
    {
	    case viergever1DOF: case diependaal1DOF: case thesis1DOF: 
        { 
	        /* basic 1DOF model */
            /* just to make sure ... there's no coupling... */
            CUTE = FALSE;
		    if ((capacitiveCircuit)||(extracellularGradient)) object[ELMO]=TRUE;
		    object[BAM]=TRUE;
		    objfft[BAM]=0;
		    channel[SM]=TRUE;
		    channel[ST]=TRUE;
		    fluidbound[BAM]=TRUE;
		    incident[BAM][SM]=TRUE;
		    incident[BAM][ST]=TRUE;
		    if (modeltype == thesis1DOF) 
            { 
			    *geo=thesisGEO;
			    *imp=thesisZ1DOF;
			    *me=noME;
		    }
		    else if (modeltype == viergever1DOF) 
            { 
			    *geo=viergeverGEO;
			    *imp=viergeverZ;
			    *me=noME;
		    }
		    else if (modeltype == diependaal1DOF) 
            { 
			    *geo=diepenGEO;
			    *imp=diepenZ;
			    *me=diepenME;
		    }
		    *cp=noCP;
		    break;
		}
	    case passive2DOF: case neely2DOF: 
        { 
	    /* basic 2DOF model */
		    /* define objects */
		    object[BAM]=TRUE;
		    object[TM]=TRUE;
		    objfft[BAM]=2;
		    objfft[TM]=1;
		    objfft[OHCSTEREO]=0;

		    /* define coupling */
		    coupled[BAM][TM]=TRUE;
		    coupled[TM][BAM]=TRUE;

		    /* define fluidboundaries and incidences */
		    channel[SM]=TRUE;
		    channel[ST]=TRUE;
		    fluidbound[BAM]=TRUE;
		    incident[BAM][SM]=TRUE;
		    incident[BAM][ST]=TRUE;
           
		    if (MULTI) {
			    fluidbound[TM]=TRUE;
			    channel[SS]=TRUE;
			    incident[TM][SM]=TRUE;
			    incident[TM][SS]=TRUE;
		    }
          
		    /* define geometry type, impedance type, middle ear function
		       and coupling type */
		    if (modeltype == passive2DOF) {
			    *imp=thesisZ2DOF;
			    *geo=thesisGEO;
			    *cp=thesis2DOFCP;
			    *me=noME;
		    }
		    else { /* neely's model */
                neelymodel = TRUE;
			    *imp=neelyZ_PSA2DOF;
			    *cp=neelyCP;
			    *geo=neelyGEO;
			    *me=neelyME;
		    }
		    break;
		}
	    case active2DOF: 
        { 
	    /* basic 2DOF model */
		    /* define objects */
		    object[BAM]=TRUE;
		    object[OHC]=TRUE;
		    object[ELMO]=FALSE;
		    objfft[BAM]=0;

		    /* define coupling */
		    coupled[BAM][OHC]=TRUE;
		    coupled[OHC][BAM]=TRUE;

		    /* define fluidboundaries and incidences */
		    channel[SM]=TRUE;
		    channel[ST]=TRUE;
		    fluidbound[BAM]=TRUE;
		    incident[BAM][SM]=TRUE;
		    incident[BAM][ST]=TRUE;
		    /* define geometry type, impedance type, middle ear function
		       and coupling type */
		    *imp=thesisZ2DOFactive;
		    *geo=thesisGEO;
		    *cp=thesis2DOFactiveCP;
		    *me=noME;
		    break;
		}
	    case active3DOF: 
        {
		    /* define objects that are integrated */
		    object[BAM]=TRUE;
		    object[OHC]=TRUE;
		    object[TM]=TRUE;
		    object[ELMO]=TRUE;

		    objfft[BAM]=0;

		    /* deiter coupling, cuticular plate coupling */
		    coupled[BAM][OHC]=TRUE;
		    coupled[OHC][BAM]=TRUE;
		    /* cilia coupling */
		    coupled[TM][BAM]=TRUE;
		    coupled[BAM][TM]=TRUE;

		    /* fluid boundaries */
		    channel[SM]=TRUE;
		    channel[ST]=TRUE;
		    fluidbound[BAM]=TRUE;
		    incident[BAM][SM]=TRUE;
		    incident[BAM][ST]=TRUE;

		    if (MULTI) {
			    channel[SS]=TRUE;
			    fluidbound[TM]=TRUE;
			    incident[TM][SM]=TRUE;
			    incident[TM][SS]=TRUE;
		    }

		    /* define geometry type, impedance type, middle ear function
		       and coupling type */
		    *geo=thesisGEO;
		    *imp=thesisZ3DOF;
		    *me=noME;
		    *cp=thesis3DOFCP;
		    break;
	    }	
	    case passive4DOF: 
        case motile4DOFcapacitive:
	    case motile4DOFextracellular: 
	    case couple4DOFcapacitive: 
	    case couple4DOFextracellular: 
	    case stiff4DOFsoma: 
	    case stiff4DOFcilia: 
	    case stiff4DOFboth: 
	    case motilestiff4DOFsoma: 
	    case motilestiff4DOFcilia: 
	    case motilestiff4DOFboth:
        case all4DOF: 
        {

		    /* define objects that are integrated */
		    object[BAM]  = TRUE;
		    object[OHC]  = TRUE;
		    object[ACRL] = TRUE;
		    object[TM]   = TRUE;

		    if ((capacitiveCircuit)||(extracellularGradient)) 
                object[ELMO]=TRUE;

		    objfft[BAM]  = 0;
		    objfft[TM]   = 1;
		    objfft[OHCSTEREO]= 2;

		    /* deiter coupling */
		    coupled[BAM][OHC]=TRUE;
		    coupled[OHC][BAM]=TRUE;

		    /* pillar foot coupling */
		    coupled[BAM][ACRL]=TRUE;
		    coupled[ACRL][BAM]=TRUE;

		    /* cilia coupling */
		    coupled[TM][ACRL]=TRUE;
		    coupled[ACRL][TM]=TRUE;

		    /* cuticular plate coupling */
		    coupled[OHC][ACRL]=TRUE;
		    coupled[ACRL][OHC]=TRUE;

		    /* adjacent coupling */
		    coupled[OHC][OHC]=TRUE;
		    coupled[ACRL][ACRL]=TRUE;

		    /* fluid boundaries */
		    channel[SM]=TRUE;
		    channel[ST]=TRUE;
		    fluidbound[BAM]=TRUE;
		    incident[BAM][SM]=TRUE;
		    incident[BAM][ST]=TRUE;

		    if (MULTI) {
			    channel[SS]=TRUE;
			    incident[BAM][SM]=FALSE;
			    incident[ACRL][SM]=TRUE;
			    fluidbound[TM]=TRUE;
			    incident[TM][SM]=TRUE;
			    incident[TM][SS]=TRUE;
		    }


		    /* define geometry type, impedance type, middle ear function
		       and coupling type */
		    *geo=thesisGEO;
		    *imp=thesisZ4DOF;
		    *me=noME;
		    *cp=thesis4DOFCP;

		    break;
	    }	
	    default: 
        { 
		    error("Watch out!  that model no longer EXISTS","topology.c");
		    break;
	    }
    }
}
/*---------------------------------------------------------------------------*/
