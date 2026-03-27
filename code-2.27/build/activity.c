#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"

/*------------------------------------------------------------------------*/
void activity(int mtype)
/*------------------------------------------------------------------------*/
{
    /* defaults */
    motileOHC               = FALSE;
    coupledOHC              = FALSE;
    capacitiveCircuit       = FALSE;
    extracellularGradient   = FALSE;
    activeCOREYstiffness    = FALSE;
    activeOHCstiffness      = FALSE;
    activeDEITERcouple      = FALSE;


    switch(mtype) {
	    case neely2DOF:  
        case diependaal1DOF:  
        case viergever1DOF: 
	    case passive2DOF: 
        case passive4DOF: 
            {
	        /* no activity */
		    break;
		    }
        case thesis1DOF: {
            capacitiveCircuit       = TRUE;
            motileOHC               = TRUE;
            coupledOHC              = FALSE; 
            break;
            }
	    case active2DOF: {
		    motileOHC               = TRUE;
		    activeOHCstiffness      = FALSE;
		    break;
		    }
	    case active3DOF: {
		    break;
		    }	
	    case motile4DOFcapacitive:  
            {
            motileOHC               = TRUE;
            coupledOHC              = TRUE;
            capacitiveCircuit       = TRUE;
		    break;
		    }
	    case motile4DOFextracellular:  
            {
            motileOHC               = TRUE;
            coupledOHC              = TRUE;
            extracellularGradient   = TRUE;
    	    break;
		    }
	    case couple4DOFcapacitive:  
            {
            motileOHC               = TRUE;
            coupledOHC              = TRUE;
            capacitiveCircuit       = TRUE;
            activeDEITERcouple      = TRUE;
		    break;
		    }
	    case couple4DOFextracellular:  
            {
            motileOHC               = TRUE;
            coupledOHC              = TRUE;
            extracellularGradient   = TRUE;
            activeDEITERcouple      = TRUE;
		    break;
		    }
        case stiff4DOFsoma:  {
            motileOHC               = TRUE;
            capacitiveCircuit       = TRUE;
            activeOHCstiffness      = TRUE;
		    break;
		    }
        case stiff4DOFcilia:  {
            activeEVANSstiffness    = TRUE;
		    break;
		    }
        case stiff4DOFboth:  {
            motileOHC               = TRUE;
            capacitiveCircuit       = TRUE;
            activeEVANSstiffness    = TRUE;
            activeOHCstiffness      = TRUE;
		    break;
		    }
        case motilestiff4DOFsoma:  {
            motileOHC               = TRUE;
            coupledOHC              = TRUE;
            capacitiveCircuit       = TRUE;
            activeOHCstiffness      = TRUE;
		    break;
		    }
        case motilestiff4DOFcilia:  {
            motileOHC               = TRUE;
            coupledOHC              = TRUE;
            capacitiveCircuit       = TRUE;
            activeEVANSstiffness    = TRUE;
		    break;
		    }
        case motilestiff4DOFboth:  {
            motileOHC               = TRUE;
            coupledOHC              = TRUE;
            capacitiveCircuit       = TRUE;
            activeEVANSstiffness    = TRUE;
            activeOHCstiffness      = TRUE;
		    break;
		    }
        case all4DOF:  {
            motileOHC               = TRUE;
            coupledOHC              = TRUE;
            capacitiveCircuit       = TRUE;
            activeEVANSstiffness    = TRUE;
            activeOHCstiffness      = TRUE;
            activeDEITERcouple      = TRUE;
		    break;
		    }

	    default: {
		    error("Activity undefined","in build/activity.c");
		    break;
	        }
    }
    if ((motileOHC)&&(!((extracellularGradient)||(capacitiveCircuit))))
        error("In activity.c,","A motile OHC is defined but OHC electrical properties are not.");
}
/*------------------------------------------------------------------------*/
