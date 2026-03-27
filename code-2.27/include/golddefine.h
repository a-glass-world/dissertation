#define DONTCARE    99
#define M_PI        3.1415926535
#define OHC_k      13.8 		/* 1.38 x 10^(-23) J/K */
#define OHC_T     273.0 	    /* temperatue in K */
#define OHC_q       1.6 		/* charge, electron */

/*------------------------------------------------------------------------*/
/*				A matrix				  */
/*------------------------------------------------------------------------*/
#define DIAG       1
#define SUBDIAG    0
#define SUPRADIAG  2

#define LDIAG      1
#define LSUBDIAG   0

#define UDIAG      0
#define USUPRADIAG 1
/*------------------------------------------------------------------------*/
/*				Files 					  */
/*------------------------------------------------------------------------*/
#define globaldir	"bin"
#define testdir		"test"
#define datadir    	"data"
#define viscousdir	"viscous"
#define clickdir	"input"

#define explainname	"/model_description"
#define snapxname	"/x"
#define snapvname	"/v"
#define snaptimename "/time"
#define fftname		"/fft"
#define fundname	"/fund"
#define scanname	"/scan"
#define measurename	"/measure"
#define debugname	"/debug"
#define ohcname		"/motion"

#define nchar		200
#define nsignal		4

#define ndebug		8
#define nfile		6	
#define nohcfile	20	

#define SNAP        0

/*------------------------------------------------------------------------*/
/*				Animal status				  */
/*------------------------------------------------------------------------*/
#define live	0
#define furo	1
#define loud	2
#define post	3

/*------------------------------------------------------------------------*/
/*				Input 					  */
/*------------------------------------------------------------------------*/
#define ntone	 2

#define TONE	 0
#define CLICK	 1
#define TTONED 	 2
#define TTONES 	 3 
#define PIP 	 4
#define CLICKFFT 5
#define NONE     6

/*------------------------------------------------------------------------*/
/*				Data					  */
/*------------------------------------------------------------------------*/
#define liveclick	    0
#define isointensity	1
#define ttoned		    2
#define ttones		    3
#define livetune	    4
#define quicklivetune   5
#define livegain	    6

#define loudclick	    7
#define loudiso		    8

#define postclick       9
#define posttune	    10

#define furoclick	    11
#define furoinout	    12
#define furoclickfft    13

#define numlclick	     4
#define numliso		    14
#define numttoned	     6
#define numttones	    10
#define numltune	    15
#define numqtune	     5
#define numgain		    13

#define numpclick	     3
#define numptune	     7

#define numfclick   	 0
#define numfio		     9
#define numfclickfft	 3

#define numloudclick	 3
#define numloudiso	    12

#define maxscan		    16 /* largest of the individual types */

/*------------------------------------------------------------------------*/
/*			Input type					  */
/*------------------------------------------------------------------------*/
#define diepInput       0 
#define neelyInput      1 
#define scanInput       2 
#define viergeverInput  3
#define userInput       4
#define fixedInput      5
/*------------------------------------------------------------------------*/
/*				General					  */
/*------------------------------------------------------------------------*/
#define FALSE	        0
#define TRUE	        1

#define clickdur        0.000782

#define length		    128
#define fftlength	    64	
#define TIMESTEPS	    512
#define nd		        2

#define stable 		    0
#define imaginary       1
#define unstable 	    2	
#define badmap	 	    3

#define iNF		        minmem	
#define minmem		    3	
#define maxfeatures	    10	

#define nobject	    	7
#define BAM		        0
#define TM		        1
#define OHC		        2
#define ACRL		    3
#define ELMO		    4	
#define OHCSTEREO       5	
#define DELL		    6

#define nchannel    	3
#define SM	        	0
#define ST	        	1
#define SS	    	    2	

#define viergever1DOF   0
#define thesis1DOF	    1
#define passive2DOF     2 
#define passive4DOF     3 
#define active2DOF      4 
#define active3DOF      5 

/* active 4DOF models in thesis */
#define motile4DOFcapacitive     6
#define motile4DOFextracellular  7
#define couple4DOFcapacitive     8
#define couple4DOFextracellular  9
#define stiff4DOFsoma           10
#define stiff4DOFcilia          11
#define stiff4DOFboth           12
#define motilestiff4DOFsoma     13
#define motilestiff4DOFcilia    14
#define motilestiff4DOFboth     15
#define all4DOF                 16

#define diependaal1DOF          20 
#define neely2DOF  	            21	

#define thesisGEO       0
#define neelyGEO        1
#define diepenGEO       2
#define chap3GEO        3
#define viergeverGEO   	4 

#define neelyME         0
#define diepenME        1 
#define noME            2 

#define thesisZ1DOF         0
#define thesisZ2DOF         1  
#define thesisZ2DOFactive   2  
#define thesisZ3DOF         3  
#define thesisZ4DOF    	    4 
#define neelyZ_PSA2DOF 	    5 
#define neelyZ_PSA1DOF      6 
#define diepenZ             7 
#define viergeverZ     	    8 

#define noCP        	    0
#define neelyCP     	    1
#define thesis2DOFCP        2 
#define thesis2DOFactiveCP  3 
#define thesis3DOFCP        4 
#define thesis4DOFCP        5 

/*------------------------------------------------------------------------*/
/*				Output					  */
/*------------------------------------------------------------------------*/
#define	DISPLACEMENT	0
#define	VELOCITY	    1	
/*------------------------------------------------------------------------*/
/*				runge kutta				  */
/*------------------------------------------------------------------------*/
#define rksteps	 	    7	
#define START	 	    0
#define CURRENT		    6	
/*------------------------------------------------------------------------*/
/*				estimates				  */
/*------------------------------------------------------------------------*/
#define RESIST          0
#define STIFF           1
#define MASS            2
/*------------------------------------------------------------------------*/
#define cACRLBAMstiff   3
#define cACRLOHCstiff   4
#define cACRLTMresist   5
#define cACRLTMstiff    6
#define cACRLACRLstiff  7
#define cOHCOHCstiff    8
#define cOHCBAMresist   9
#define cOHCBAMstiff   10

#define ncoef		    2
#define nparm		    11
#define nestimates	    nparm*ncoef*16	
