

/*------------------------------------------------------------------------*/
/* 				GLOBALS					  */
/*------------------------------------------------------------------------*/
/* 			upper-level control variables 			  */
/*------------------------------------------------------------------------*/
int CUTE; 		/* toggle to include coupling */
int animal_status; 	/* describe animal state for input constraints */
int MULTI;		/* define strategy used for fluid derivatives */
int viscosity;  /* if TRUE: strategy is for viscous solutions */
int vary; 		/* use the variable stepper for the Runge-Kutta */
int estimateactive;     /* estimate the active parameters vs passive */
int constantactive;     /* constant active parameters */
int RECIO;		/* determines if click input is from file or pure */
extern int CUTE, animal_status, RECIO, constantactive, estimateactive;
extern int MULTI, viscosity, vary;


/*------------------------------------------------------------------------*/
/* 				indices 				  */
/*------------------------------------------------------------------------*/
int nmodel; 		/* current model */
int nscanmodel;		/* current input paradigm */
extern int nmodel, nscanmodel;


/*------------------------------------------------------------------------*/
/* 				objects 	 			  */
/*------------------------------------------------------------------------*/
int actmem; 		/* number of objects using memory during analysis */
int objfft[nobject]; 	/* objects analyzed using FFT */
int object[nobject];	/* objects in model */
extern int actmem, objfft[nobject], object[nobject];


/*------------------------------------------------------------------------*/
/* 			estimate control variables 			  */
/*------------------------------------------------------------------------*/
int pauses; 		/* pause at various points in estimation strategy */
int newparm; 		/* if TRUE: new parameter type has been defined */
int estimating; 	/* if TRUE: estimate parameters; else use constant set */
int calibrating;	/* if TRUE: first pass, observation index unknown */
extern int pauses, newparm, estimating, calibrating; 


/*------------------------------------------------------------------------*/
/*			estimation structs  				  */
/*------------------------------------------------------------------------*/
struct scanstruct response[maxscan][nestimates][nobject]; /* model responses */
struct scanstruct success[maxscan][nestimates]; /* model success estimate */
extern struct scanstruct response[maxscan][nestimates][nobject], 
success[maxscan][nestimates]; /* model success estimate */


/*------------------------------------------------------------------------*/
/* 			estimation constraints 		 		  */
/*------------------------------------------------------------------------*/
struct scanstruct data[maxscan];/* data constraints */
extern struct scanstruct data[maxscan];

double lomap, himap; 		/* frequency place map restraints */
double freqdiff;		    /* frequency place map success estimate */
extern double lomap, himap, freqdiff;


/*------------------------------------------------------------------------*/
/* 				input variables 			  */
/*------------------------------------------------------------------------*/
double freq[ntone]; 	/* input tone array (a spectra) */
extern double freq[ntone];

int intone; 			/* number of tones in the input spectra */
int scanning; 			/* if TRUE: use range of input paradigms */
int ranging; 			/* if TRUE: use range of model parameters */
int inputtype;			/* type of input used (if not scanning) */
extern int intone, scanning;


/*------------------------------------------------------------------------*/
/* 				parameters 				  */
/*------------------------------------------------------------------------*/
double cl;	 		        /* the length of the cochlea in mm. */
double mass[nobject][2];	/* mass of object where the first term in */
				            /* the array is the constant term and the */
				            /* second is the exponential. */
double resist[nobject][2];	/* resistance of object */
double stiff[nobject][2];	/* stiffness of object */
extern double cl, mass[nobject][2], resist[nobject][2], stiff[nobject][2];	


/*------------------------------------------------------------------------*/
/* 			coupling parameters 				  */
/*------------------------------------------------------------------------*/
double cmass[nobject][nobject][2]; 	    /* mass coupling */
double cresist[nobject][nobject][2]; 	/* resistive coupling */
double cstiff[nobject][nobject][2]; 	/* stiffness coupling */
extern double cmass[nobject][nobject][2], cresist[nobject][nobject][2],
cstiff[nobject][nobject][2];


/*------------------------------------------------------------------------*/
/* 				output  				  */
/*------------------------------------------------------------------------*/
char outputdir[nchar]; 		/* top directory where data files wind up */
char currentdir[nchar]; 	/* default name for data files */
extern char outputdir[nchar], currentdir[nchar]; 		

struct fileinfo ohcfile[nohcfile]; /* these are the ohc motility files */
struct fileinfo explain;	/* file describing strategy and parameters */
struct fileinfo estimatefile;	/* estimate convergence, global response */
				                /* features, success measures, stability */
struct fileinfo rrr; /* hey, its a mystery! BRRRRRRR!!!*/
/* no actually its just a quick peek at the parameter estimation in the
   current directory - gives the features and their success per parameter */
extern struct fileinfo explain,estimatefile,rrr,ohcfile[nohcfile];

/* snapx (2d = time * displacement) */
struct fileinfo snapx[nobject],snapv[nobject];
/* clickx (2d = loc * displacement) */
struct fileinfo clickx[nobject],clickv[nobject], clickt[nobject];
extern struct fileinfo snapx[nobject],snapv[nobject],clickx[nobject],clickv[nobject], clickt[nobject];

double snapt; 			/* snap time - start time for time-domain */
extern double snapt;	/* picture of the state-space response */


/*------------------------------------------------------------------------*/
/* 				fluid equations			  */
/*------------------------------------------------------------------------*/
struct motion P[rksteps];
extern struct motion P[rksteps];

struct motion Y[rksteps][nd][nobject], D[rksteps][nd][nobject];
extern struct motion Y[rksteps][nd][nobject], D[rksteps][nd][nobject];

struct motion activeC[rksteps][nd][nobject];
extern struct motion activeC[rksteps][nd][nobject];

int YF, DF, iO;
extern int YF, DF, iO;

double mesh[length+1];
extern double mesh[length+1];

double finput, Gc[length];
extern double finput, Gc[length];

double Gme,totsignal,Gwindow;
extern double Gme,totsignal,Gwindow;

double G[length], K[length], PHI[length];
extern double G[length], K[length], PHI[length];

int obsindex;
extern obsindex;

double observe;
extern double observe;

int fftstarted;
extern int fftstarted;

int settype;
extern int settype;

double ohcscale;
extern double ohcscale;

double hifreq, lofreq;
extern double hifreq, lofreq;

double weight[maxfeatures];
extern double weight[maxfeatures];

double S[nobject][length], M[nobject][length], R[nobject][length];
extern double S[nobject][length], M[nobject][length], R[nobject][length];

double S0[length]; 
extern double S0[length];

double cS[nobject][nobject][length], cR[nobject][nobject][length];
extern double cS[nobject][nobject][length], cR[nobject][nobject][length];

double damp[nobject][length], lambda[nobject][length];
extern double damp[nobject][length], lambda[nobject][length];

double maxld[length], minld[length], rangeld[length];
extern double maxld[length], minld[length], rangeld[length];

struct perturb sensitive;
extern struct perturb sensitive;

int fluidbound[nobject];
extern int fluidbound[nobject];

double L[nobject][length][2], U[nobject][length][2];
extern double L[nobject][length][2], U[nobject][length][2];

double incident[nobject][nchannel],density, area[nchannel][length];
extern double incident[nobject][nchannel],density,area[nchannel][length];

double delta[length+1], vgamma[nobject][length], alpha[nobject][length];
extern double delta[length+1],vgamma[nobject][length], alpha[nobject][length];

double gain[length];
extern double gain[length];

int activeEVANSstiffness, activeCOREYstiffness, activeOHCstiffness, activeDEITERcouple;
extern int activeEVANSstiffness, activeCOREYstiffness, activeOHCstiffness, activeDEITERcouple;

int motileOHC, coupledOHC, capacitiveCircuit, extracellularGradient;
extern int motileOHC, coupledOHC, capacitiveCircuit, extracellularGradient;


/**************************************************************************/
/*                              ohc constants                             */
/**************************************************************************/
/* for resistive and capactive membrane */
double Vb0[length], V0, Vm, Ca0, Ca1, Cb0, Cb1, Rb0, Rb1;
extern double Vb0[length], V0, Vm, Ca0, Ca1, Cb0, Cb1, Rb0, Rb1;

/* stereocilia stiffness */
double sp0,sp1,sa,sl;
extern double sp0,sp1,sa,sl;

/* channel conductance */
double gc,gleak,stereo_N0,stereo_N1,stereo_N[length];
extern double gc,gleak,stereo_N0,stereo_N1,stereo_N[length];

/* probabilities */
double zP, X0;
extern double zP, X0;

double keepzP;
extern double keepzP;

/* ohc geometry */
double ohc_r;
extern double ohc_r;
double stereo_h0, stereo_h1, stereo_h[length];
extern double stereo_h0, stereo_h1, stereo_h[length];
double ohc_l0, ohc_l1, ohc_l[length];
extern double ohc_l0, ohc_l1, ohc_l[length];
double apical_a0, apical_a1, apical_a[length];
extern double apical_a0, apical_a1, apical_a[length];
double basolateral_a0, basolateral_a1, basolateral_a[length];
extern double basolateral_a0, basolateral_a1, basolateral_a[length];


/**************************************************************************/

double T2mRm,Tm2mm_ms, Tm, T2m, Am, Rm, Sm;
extern double T2mRm,Tm2mm_ms, Tm, T2m, Am, Rm, Sm;

int channel[nchannel], coupled[nobject][nobject], MIDEAR;
extern int channel[nchannel], coupled[nobject][nobject], MIDEAR;

int gaindefined;
extern int gaindefined;

double gain_a11, gain_a12;
extern double gain_a11, gain_a12;

int  signaltype;
extern int  signaltype;

double level[ntone], dB[ntone], scale;
extern double level[ntone], dB[ntone], scale;

int meshsize;
extern int meshsize;

double width[nobject][2],height[nchannel][3];
extern double width[nobject][2],height[nchannel][3];

double ZweigN;
extern double ZweigN;

struct fileinfo fftfile[nfile][minmem];
extern struct fileinfo fftfile[nfile][minmem];

struct fileinfo fundfile[nfile][minmem];
extern struct fileinfo fundfile[nfile][minmem];

struct fileinfo debug[nobject][ndebug];
extern struct fileinfo debug[nobject][ndebug];

struct fileinfo infreqfile, indBfile;
extern struct fileinfo infreqfile, indBfile;

char currentdirname[nchar];
extern char currentdirname[nchar];

/* used in buildestimate */

double bestparm[nparm][ncoef][nobject+1];
extern double bestparm[nparm][ncoef][nobject+1];

int plate[nobject], UC[nobject], LC[nobject];
extern int plate[nobject], UC[nobject], LC[nobject];

int cycle;
extern int cycle;

int sense_cycle[nparm][ncoef];
extern int sense_cycle[nparm][ncoef];

double currentparm;
extern double currentparm;

double range[2][nparm][nobject+1][2];
extern double range[2][nparm][nobject+1][2];

double keepstiff[nobject][2], keepmass[nobject][2], keepresist[nobject][2]; 
extern double keepstiff[nobject][2],keepmass[nobject][2],keepresist[nobject][2]; 

double keepcstiff[nobject][nobject][2], keepcresist[nobject][nobject][2];
extern double keepcstiff[nobject][nobject][2], keepcresist[nobject][nobject][2];

struct td tds;
extern struct td tds;

int datatype;
extern int datatype;

struct fft_in ffIn[iNF][fftlength];
extern struct fft_in ffIn[iNF][fftlength];

struct fft_state fftFreq[iNF][fftlength];
extern struct fft_state fftFreq[iNF][fftlength];

double Gamma0[nchannel][length], Gamma1[nchannel][length] ;
extern double Gamma0[nchannel][length], Gamma1[nchannel][length] ;

double beta[nobject][length];
extern double beta[nobject][length];

struct vstepper vs;
extern struct vstepper vs;

int peakindex;
extern int peakindex;

double shape_factor_a[length];
extern double shape_factor_a[length];

double restingRa[length];
extern double restingRa[length];

double FFT_max_mag[iNF];
extern double FFT_max_mag[iNF];

int FFT_max_comp[iNF] ,FFT_max_loc[iNF];
extern int FFT_max_comp[iNF] ,FFT_max_loc[iNF];

int neelymodel;
extern int neelymodel;

int estimator_bound[6];
extern int estimator_bound[6];

double staticS[length];
extern double staticS[length];

/* new net variables to take care of the slow exp */
double meshRb[length], meshCa0[length], meshCb0[length];
extern double meshRb[length], meshCa0[length], meshCb0[length];
double expb0;
extern double expb0;