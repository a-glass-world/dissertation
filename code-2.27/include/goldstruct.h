/*-----------------------------------------------------------------------*/
/*			goldstruct.h					 */
/*-----------------------------------------------------------------------*/
struct fft_state {
	double freq[TIMESTEPS];
	double real[TIMESTEPS];
	double im[TIMESTEPS];
	double mag[TIMESTEPS];
	double phase[TIMESTEPS];
	};
/*-----------------------------------------------------------------------*/
struct stacktype {
	double x;
	double p;
	}; 
/*-----------------------------------------------------------------------*/
struct fft_in {
	double time[TIMESTEPS];
	double real[TIMESTEPS];
	double im[TIMESTEPS];
	}; 
/*-----------------------------------------------------------------------*/
struct fileinfo {
	int open;
	FILE *fp;
	char c[nchar]; 
	}; 
/*-----------------------------------------------------------------------*/
struct scanstruct {
	double observe;
	int stability;
	int obsindex;
	int maxindex;
	int fundindex;
	int modelnum;
	int inputnum;
	int signal;
	int dattype;
	int nfeature;
	int cal;
	int newpm;
	int estcoef;
	int estparm;
	int estobj;
	double est;
	double scan_norm_tot;
	double norm;
	double freq;
	double dB;
	double feature[maxfeatures];
	double diff[maxfeatures];
	double deriv[maxfeatures];
	double intval;
	} ;
/*-------------------------------------------------------------------------*/
struct perturb {
	double nmod[nparm][ncoef];
	double featinp[nparm][ncoef][maxscan];
	double feature[nparm][ncoef][maxfeatures];
	} ;
/*-------------------------------------------------------------------------*/
struct td {
	int runfft;
	int outputType;
	int print3Dfft;
	int iMAXSTEP;
	int defined;
	int iskip;
	int istable;
	int dps;
	double DOF;
	double starttime;
	double endtime;
	double timestep;
	double normalstep;
	double fft_timestep;
	double cycletime;
	double fftstepspercycle;
	double fftcyclesperstep;
	double rkstepspercycle;
	double rkcyclesperstep;
	double stabletime;
	double minstep;
	};
/*-------------------------------------------------------------------------*/
struct vstepper {
	double maxerr;
	double minerr;
	double decrement;
	int startover;
	};
/*-------------------------------------------------------------------------*/
struct motion {
	double stapes;
	double window;
	double loc[length];
	};
/*-------------------------------------------------------------------------*/
struct corey_complex {
	double re;
	double im;
	double ma;
	double ph;
	};
/*-------------------------------------------------------------------------*/
struct topology {
	double x;
	};
/*-------------------------------------------------------------------------*/
struct gradient {
	double diff[maxfeatures];
	double sens[maxfeatures];
	};
/*-------------------------------------------------------------------------*/
