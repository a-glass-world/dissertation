
/* in planprint.c */
void list_topology();
void list_geometry();
void list_impedance();
void list_coupling();
void list_middleear();
void list_input();
void list_ohc();
void nameit(
            int o, 
            FILE *fp); 
void namechannel();

/* utility functions */
void error(char message1[200],char message2[200]);

int query();
void shortsummary();


/* in A.c */
void buildA(int DEBUG);
void constructA(double A[length][3], int obj);
void buildLU(double A[length][3], int obj);
void debugA(double A[length][3], int obj);

/* in activity.c */
void activity(int mtype);

/* in buildestimate.c */
int buildestimate(
                  int stability, 
                  int mtype, 
                  int *px, 
                  int *c, 
                  int *o);

void getbestestimates();

void getbestestimate(
                     int p, 
                     int c, 
                     int o);

void savebestestimate(
                      int p, 
                      int c, 
                      int o);

void getcoupledobjs(
                      int p, 
                      int *co1, 
                      int *co2);

void keepparms();

int changeparameter(
                    int mtype, 
                    int *px,
                    int *c,
                    int *o,
                    int *newcycle);

int changeestimate(
                   int ptype,
                   int pcoef,
                   int obj,
                   int st);

double estimator(
                 int ptype, 
                 int obj,
                 int pcoef, 
                 int newstack, 
                 int stability, 
                 int *done);

int sort(
         struct stacktype s[nestimates],
         int n);

double GetPlateMaterial(
                        int    ptype,
                        int    pcoef,
                        double parameter,
                        int    obj);

double SetPlateMaterial(
                        int    ptype,
                        int    pcoef,
                        double inparm,
                        int    obj);
void setplateranges();

void printMOREstuff(
                    FILE *fp,
                    int ptype,
                    int pcoef, 
                    int obj);

double ui_returnvalue(
                    int ptype,
                    int pcoef, 
                    int obj);

char *ui_returnname(int ptype);

double resetparm(
                 int ptype,
                 int pcoef,
                 int obj);

int ui_estimate(
                int p,
                int c,
                int o,
                int ui_type);
void getbounds(
        int *mino,
        int *maxo,
        int *minp,
        int *maxp,
        int *minc,
        int *maxc);

void setbounds(
        int mino,
        int maxo,
        int minp,
        int maxp,
        int minc,
        int maxc);


/* in buildohc.c */
void initOHC();
void OHCgeometry();

/* in buildresponse.c */
void comment(int n);

void buildresponse(
                   int n,
                   int p,
                   int c,
                   int o,
                   int np,
                   int stability,
                   int obsindex);
void buildclickfftfeat(
                       int n,
                       int m,
                       int i,
                       double calfreq,
                       int obsindex);
void buildclickfeat(
                    int n,
                    int m,
                    int i,
                    int fundindex,
                    int obsindex);
void buildtimefeat(
                   int i,
                   double *onset,
                   double *decay,
                   double *duration,
                   double *peakvel,
                   double *numcyles,
                   int obsindex);
void buildresponsetop(
                      int n,
                      int m,
                      int i,
                      int fundindex,
                      int obsindex,
                      int maxindex,
                      int ptype,
                      int obj,
                      int coef,
                      int np,
                      int stability);
void buildisofeat(
                  int n,
                  int m,
                  int i,
                  int fundindex,
                  int obsindex);

void buildtunefeat(
                   int n,
                   int m,
                   int i,
                   int fundindex,
                   int obsindex);
void buildtunefeat(
                   int n,
                   int m,
                   int i,
                   int fundindex,
                   int obsindex);
void buildgainfeat(
                   int n,
                   int m,
                   int i,
                   int fundindex,
                   int obsindex);
double measuregain(
                   double velocity,
                   double dBSPL);
void getfourierindex(
                     int *index,
                     int obj,
                     double comp);

void getplaceindex(
                   int k,
                   int *index,
                   int obj);



/* in buildsuccess.c */
void buildsuccess(
                  int ninput,
                  int dtype,
                  int stability);

int makesuccessfeature(
                        int n,
                        int m);
void resetsensitivity(
                      int p,
                      int c,
                      int n,
                      int m);
int calculable(
               int n,
               int m);

int sensitivity_calculated(
                           int p,
                           int c);
void setsensitive(
                  int n,
                  int m);

void constructweights(
                      int n,
                      int m,
                      int *nweight);

void measuresuccess(
                    int n,
                    int m,
                    int dtype,
                    int nweight);

void measuretrajectory(int m);

void hyperplane(struct topology intersect);

void hessian(
            struct topology intersect,
            struct topology partitions);

void bifurcate(struct topology partitions);


/* in componentprint.c */
void printFund(int   m);

void printfundvalue(
                    int i,
                    int j,
                    int k);

int makefundname(
                 int   i,
                 int   j,
                 int   m);

void getobjectname(
             char objname[20],
             int obj);

void getextension(
             char ext[5],
             int exttype);


/* in coupling.c */
void coupling(int cp);

/* in dell.c */
void dell(
          double t,
          int YIN,
          int index);
void setzerocouple(
                   int YIN,
                   int n);
void setcouple(
               int YIN,
               double dL,
               double t,
               int index,
               int n);
void changecoupling(
                    double x,
                    double t,
                    int i,
                    int index);
double ohcshape(
                double ptn,
                int i);
void setOHCstiffness(
                     double dL,
                     int i);

/* in elmo.c */
void elmo(
          double t,
          int YIN, 
          int DOUT, 
          int index);
double deflect(
               double t,
               int i,
               int YIN);
double apical_channel_probability(double x);
double bundle_stiffness(
                        int YIN,
                        double pa,
                        int i,
                        double x);
double apical_resistance(
                         double pa,
                         int i);

/* in error.c */
void closeoldmodel();

/* in fft.c */
void fft(
         int iCount,
         int *bDone,
         double inTime);

void ConstructFFTinput(
                       int iCount,
                       int YIN,
                       double inTime);
void runFFT(double roll[TIMESTEPS][minmem]);
void setdata(
             int j,
             int x,
             double dat[TIMESTEPS*2 +1]);
void FFTguts(
             double dat[2*TIMESTEPS+1],
             int nn, 
             int isign);
void setffFreq(
               int j,
               int x,
               double t[TIMESTEPS*2 +1],
               double roll[TIMESTEPS][minmem]);

/* in printfft.c */
void printFFT(int m);
void printvalue(
                int i,
                int j,
                int k);
void gnuprintvalue(
                int i,
                int j,
                int k,
                int fundindex);
int makename(
             int i,
             int j,
             int m);



/* in fitFP.c */
void fitFP(
           double low,
           double high,
           double apex,
           int obj);
double fpx(
           double s,
           double m,
           double r,
           int *stabi);

/* in fluiddebug.c */
void fluiddebug(
                double time,
                int x, 
                int y,
                int OB);

/* in fluids.c */
int fluids(
           int i,
           double time,
           int in,
           int out,
           int DEBUG);
void GetSignal(double time);
void CalculateG();
double getactivecouple(
                 int i,
                 int comp);
double getpassivecouple(
                 int i,
                 int comp);
void CalculateGc();
double BoundaryConditions();
void QuickDerive();
void CalculateK(double insignal);
void CalculateP();
void DeriveStapes();
void Derive();
int CheckUnstable();

/* in furoinput.c */
void furoscan(
              int n,
              int measuretype);
void furoclickinput(int n);
void furocfftinput(int n);
void furoioinput(int n);

/* in geometry.c */
void geometry(int GEOtype);

/* in impedance.c */
void impedance(int Ztype);
void CPimpedance();
void BAMimpedance();
void TMimpedance();
void ACRLimpedance();
void OHCimpedance();
void ELMOimpedance();

/* in initstrategy.c */
void initTD();
void resetTD();	

/* in inputanddata.c */
void defineinput(
                 int n, 
                 int inptype,
                 int datatype);
void calibrate(
               double obs,
               int n);


/* in liveinput.c */
void livescan(
              int n,
              int measuretype);
void liveclickinput(int n);
void liveisoinput(int n);
void livettonedinput(int n);
void livettonesinput(int n);
void quicktune(int n);
void tuningcurve(int n);
void gaininput(int n);

/* in loudinput.c */
void loudscan(
              int n,
              int measuretype);
void loudisoinput(int n);
void loudclickinput(int n);

/* in main.c */
void initUpperControls();
void mainEventLoop();
void systeminitplan(int modelnumber);
void userinitplan(
                  int *datatype,
                  int *modelnumber);

void appendDname(
                 char putdir[nchar],
                 int dtype,
                 int addpunc,
                 char getdir[nchar]);

int model(
          int starting,
          int obs,
          int modeltype,
          int nm,
          int coef);

int decide(
           int modeltype, 
           int *ninput, 
           int edone,
           int stability,
           int p,
           int c,
           int o);

void openestimate(int dtype);

void sortresponses(
                   int totmod,
                   int datatype);
void constraints(
                 int *stability,
                 int obs);
void printpfile(int modeltype);
void printparm(FILE *fp,
               int k);
void printcparm(
                FILE *fp,
                int obj,
                int cobj);
void printOHCparm(FILE *fp);
int setobs(double infreq, int mtype);

/* in math.c */
double exp10(double power);

struct corey_complex cmult(
                           struct corey_complex a, 
                           struct corey_complex b);
struct corey_complex cdiv(
                          struct corey_complex a, 
                          struct corey_complex b);

/* in mfluids.c */
int Mfluids(
            int i,
            double time,
            int in,
            int out,
            int DEBUG);
double MBoundaryConditions(); 
double calcBC(int chan);
int MCalculatePHI();
void doPHI(
           int chan,
           int obj);
void MCalculateK(double insignal);
void MDeriveRWindow();


/* in middleear.c */
void middleear(int MEtype);

/* in net.c */
int buildnet(
             int modeltype,
             int starting);
void buildmesh();
void buildgeometry();
int buildimpedance();
void buildcoupling();
void buildalpha(int modeltype);
void buildmisc(); 

/* in netprint.c */
void netprint(FILE *fp);       
void printcZ(
             FILE *fp,
             double s1,
             double s2,
             double r1,
             double r2);

/* in ohcadmin.c */
void tellOHCparms();
void writeOHCfile(
                  double x,
                  double y,
                  char n[],
                  int fn);

void closeOHCfiles();
void openOHCfiles(
                  char n[],
                  int i);
void makeohcname(
                  char n[],
                  int i);

/* openoutput.c */
void openoutput(
                int n,
                int m,
                int c,
                int nm,
                int aut);
void makemodeldir(
                  int n,
                  int nm, 
                  int c,
                  int rindex,
                  char modeldirname[nchar],
                  int aut);
void closeoldmodel();

/* planprint.c */
void planprint(
               FILE *fp,
               int n,
               int mtype);
void list_topology(FILE *fp);
void list_geometry(FILE *fp);
void list_impedance(FILE *fp);
void list_middleear(FILE *fp);
void list_coupling(FILE *fp);
void list_input(
                FILE *fp,
                int n,
                int m);
void list_ohc(FILE *fp);
void nameit(
            int i,
            FILE *fp);
void namechannel(
                 int i,
                 FILE *fp);

/* postinput.c */
void postscan(
              int n,
              int measuretype);
void posttuningcurve(int n);
void postclickinput(int n);

/* query.c */
int query(
          char cMsg[150], 
          int iLower, 
          int iUpper);

/* readclick.c */
double getclick(double ntime);
void openclick();
int readclick(FILE *fp);

/* rk4.c */
int rk4(
        int j,
        double time,
        double step);
int derive(
           int obj,
           double time,
           double step,
           int IN,
           int OUT,
           int indexstep);
void varystep(
              struct motion x,
              struct motion y,
              double time);
void varyTD(int varytype);

/* scanprint.c */
void printscan(
               int n,
               int obj);
FILE *openscanfiles(
                    int j,
                    int m);
void showresponse(
                  int i,
                  int j,
                  int obj,
                  FILE *fp);
void makescanname(
                  int obj,
                  int m);

/* snapprint.c */
void snapprint(
               int YIN,
               double time,
               int inloc,
               int st,
               int nm, 
               int newx);
struct fileinfo opensnapfiles(
                              int nm,
                              int obj,
                              char it[10]);

/* solve.c */
int solve(
          int newoutput,
          int nm);
void setoutput(
               int count,
               double time,
               int newx,
               int nm);
int stopTD(
           double time,
           int *done);
void store(
           int IN,
           int OUT);
void initConditions();
void initZeroConditions(int i);
void initELMOConditions();


/* strategyprint.c */
void strategyprint(FILE *fp);

/* summary.c */
void quickfile(
               int dtype,
               int inp,
               int mod);
void newfile(
             int num,
             char text1[],
             char text2[],
             int showmodelnumber);
void summary(
             FILE *fp,
             int dtype,
             int inp,
             int mod,
             int solved);
void summaryscan(
                 FILE *pf,
                 struct scanstruct x,
                 int mod,
                 int dtype,
                 int flag);
void shortsummary(
                  FILE *fp,
                  int obj,
                  int parmtype);

/* topology.c */
void topology(
              int modeltype,
              int *geo,
              int *imp,
              int *me,
              int *cp);

/* trans.c */
void resistivemembrane(
                     double Ra,
                     double x,
                     int YIN,
                     int DOUT,
                     double t,
                     int i);

void capacitivemembrane(
                     double Ra,
                     double x,
                     int YIN,
                     int DOUT,
                     double t,
                     int i);

double VbRest(
              double Ra,
              double Rb);

/* viscous.c */
void buildGamma(int chan);
struct corey_complex readviscoustable(
                                      FILE *fp,
                                      double infreq);
void initvgamma(int obj);
void buildvgamma(int obj);
void calcwnumber(
                 int obj,
                 double wn[length]);
void checkgamma(
                int YF,
                int iO,
                double time);

