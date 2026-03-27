/**************************************************************************
 *  buildestimate.c
 **************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gold.h"		

#define debugging FALSE

#define NOTFOUND	0	
#define FOUND		1
#define lo		0
#define hi		1


/*------------------------------------------------------------------------*/
int buildestimate(
                  int stability, 
                  int mtype, 
                  int *p, 
                  int *c, 
                  int *o)
/*------------------------------------------------------------------------*/
{
int done = FALSE;
int minp, maxp, minc, maxc, mino, maxo;
static int newcycle;

/* first model - initialize  ranges, indices and parameters */
if (nmodel == 0) 
{
    mino = BAM+1;
    maxo = nobject;
    minp = 3;
    maxp = 10;
    minc = 0;
    maxc = 1;
    setbounds(mino,maxo,minp,maxp,minc,maxc);

	/* get the indices for the first parameter to be estimated */
	done = changeparameter(mtype, p, c, o, &newcycle);
    if ((*p != minp) || (*c != minc) || (*o != mino))
        error("This is a stupid mistake.","in buildestimate()");

	if (done) return(TRUE);

	/* at the very beginning, set up the ranges for plate equations */
	setplateranges(mtype);

	/* save the initial parameter set */
	keepparms();
}

/* estimate the parameter */
newparm = changeestimate(*p, *c, *o, stability);

/* when the estimator is done with this parameter, changeparm returns 
   newparm = TRUE and sets the current parameter in the model to the 
   best value.  This is what we want to save - but we don't want to
   solve the model with this value ...*/
if (newparm)
	/* save the last estimated value for the parameter */
	savebestestimate(*p, *c, *o);

while (newparm) 
{
	/* then reset parameter to its initial value for the cycle */
	(void) resetparm(*p, *c, *o);

	/* get a new parameter - because we don't want to solve the
	   old one anymore ... */
	done = changeparameter(mtype, p, c, o, &newcycle);
	if (done) 
	{
		return(TRUE);
	}
	else
	{
	   if (newcycle) 
	   {
		getbestestimates();
		keepparms(); 	/* save the best estimates */
	   }
	   /* run the estimator again so we get a new parameter */
	   newparm = changeestimate(*p, *c, *o, stability);
	}
}

return(done);
}
/*------------------------------------------------------------------------*/
void getbestestimates()
/*------------------------------------------------------------------------*/
{
    int p, c, o;
    int minp, maxp, mino, maxo, minc, maxc;

    getbounds(&mino, &maxo, &minp, &maxp, &minc, &maxc);

    for (o = mino; o < maxo; o++) 
    { 
        if (object[o]) 
        {
            for (p = minp; p <= maxp; p++) 
            {
	            for (c = minc; c <= maxc; c++) 
                {
	            getbestestimate(p, c, o);
                }
            }
        }
    }
}
/*------------------------------------------------------------------------*/
void getbestestimate(
                     int p, 
                     int c, 
                     int o)
/*------------------------------------------------------------------------*/
{
    int co1, co2;
    getcoupledobjs(p,&co1,&co2);

    switch(p) 
    {
	    case MASS:      {mass[o][c]    = bestparm[p][c][o]; break;}
	    case STIFF:     {stiff[o][c]   = bestparm[p][c][o]; break;}
	    case RESIST:    {resist[o][c]  = bestparm[p][c][o]; break;}

        /* coupled parameters */
        case cACRLTMstiff:  
        case cACRLOHCstiff: 
        case cACRLBAMstiff: 
        case cACRLACRLstiff:
        case cOHCOHCstiff:  
        case cOHCBAMstiff:  {cstiff[co1][co2][c]  = cstiff[co2][co1][c]  = bestparm[p][c][nobject]; break;}

        case cACRLTMresist: 
        case cOHCBAMresist: {cresist[co1][co2][c] = cresist[co2][co1][c] = bestparm[p][c][nobject]; break;}

        default: 
        { 
		    error("Parameter type undefined.", "in savebestestimate()");  
        }
    }

}
/*------------------------------------------------------------------------*/
void savebestestimate(
                      int p, 
                      int c, 
                      int o)
/*------------------------------------------------------------------------*/
{
    int co1, co2;
    getcoupledobjs(p,&co1,&co2);

    switch(p) 
    {
	    case MASS:          {bestparm[p][c][o] = mass[o][c];   break;}
	    case STIFF:         {bestparm[p][c][o] = stiff[o][c];  break;}
	    case RESIST:        {bestparm[p][c][o] = resist[o][c]; break;}

        case cACRLTMstiff:  
        case cACRLOHCstiff: 
        case cACRLBAMstiff: 
        case cACRLACRLstiff:
        case cOHCOHCstiff:  
        case cOHCBAMstiff:  {bestparm[p][c][nobject] = cstiff[co1][co2][c];  break;}

        case cACRLTMresist: 
        case cOHCBAMresist: {bestparm[p][c][nobject] = cresist[co1][co2][c]; break;}

	    default: 
        { 
		    error("Parameter type undefined.", "in savebestestimate()");  
        }
    }

}
/*------------------------------------------------------------------------*/
void getcoupledobjs(
                      int p, 
                      int *co1, 
                      int *co2)
/*------------------------------------------------------------------------*/
{
    switch(p)
    {
        case cACRLTMstiff:  {*co1 = ACRL; *co2 = TM;   break;}
        case cACRLOHCstiff: {*co1 = ACRL; *co2 = OHC;  break;}
        case cACRLBAMstiff: {*co1 = ACRL; *co2 = BAM;  break;}
        case cACRLACRLstiff:{*co1 = ACRL; *co2 = ACRL; break;}
        case cOHCOHCstiff:  {*co1 = OHC;  *co2 = OHC;  break;}
        case cOHCBAMstiff:  {*co1 = OHC;  *co2 = BAM;  break;}
        case cACRLTMresist: {*co1 = ACRL; *co2 = TM;   break;}
        case cOHCBAMresist: {*co1 = OHC;  *co2 = BAM;  break;}

	    default: {break;}
    }

}	
/*------------------------------------------------------------------------*/
void keepparms()
/*------------------------------------------------------------------------*/
{
    int i,j,k;

	for (i=0;i<nobject;i++) 
    {
	    for (j=0;j<2;j++) 
        {
		    keepmass[i][j]  = mass[i][j];
		    keepresist[i][j]= resist[i][j];
		    keepstiff[i][j] = stiff[i][j];
		}
	}

	for (i=0;i<nobject;i++) 
    {
	    for (j=0;j<nobject;j++) 
        {
		    if (coupled[i][j]) 
            {
	    	    for (k=0;k<2;k++) 
                {
					keepcresist[i][j][k] = cresist[i][j][k];
					keepcstiff[i][j][k]  = cstiff[i][j][k];
				} 
			}
		}
    }
}
/*------------------------------------------------------------------------*/
double resetparm(
                 int ptype,
                 int pcoef,
                 int obj)
/*------------------------------------------------------------------------*/
{
    int co1, co2;

    getcoupledobjs(ptype,&co1,&co2);
    switch(ptype) 
    {
	    case MASS:  
        {
            mass[obj][pcoef] = keepmass[obj][pcoef];
            return(SetPlateMaterial(ptype,pcoef,mass[obj][pcoef],obj));
        } 
	    case STIFF:  
        {
            stiff[obj][pcoef] = keepstiff[obj][pcoef];
            return(SetPlateMaterial(ptype,pcoef,stiff[obj][pcoef],obj));
        } 
	    case RESIST:  
        {
            resist[obj][pcoef] = keepresist[obj][pcoef];
            return(SetPlateMaterial(ptype,pcoef,resist[obj][pcoef],obj));
        } 
        case cACRLTMstiff:  
        case cACRLOHCstiff: 
        case cACRLBAMstiff: 
        case cACRLACRLstiff:
        case cOHCOHCstiff:  
        case cOHCBAMstiff:
        {
            cstiff[co1][co2][pcoef] = cstiff[co2][co1][pcoef] = keepcstiff[co1][co2][pcoef];
            return(cstiff[co1][co2][pcoef]);
        }
        case cACRLTMresist: 
        case cOHCBAMresist: 
        {
            cresist[co1][co2][pcoef] = cresist[co2][co1][pcoef] = keepcresist[co1][co2][pcoef];
            return(cresist[co1][co2][pcoef]);
        }
	    default: 
        { 
            error("In buildestimate.c, in resetparm()","parameter type is undefined.");
	    }
    }
}
/*------------------------------------------------------------------------*/
int changeparameter(
                    int mtype, 
                    int *p,
                    int *c,
                    int *o,
                    int *newcycle)
/*------------------------------------------------------------------------*/
{
    static int started, matparmdone;
    int getout;
    int quickcycle = FALSE; 
    int getobj = *o;
    int minp, maxp, mino, maxo, minc, maxc;

    getbounds(&mino, &maxo, &minp, &maxp, &minc, &maxc);

	if (quickcycle) 
	    maxp = 1;
	
	*newcycle = FALSE;

	if (!started) 
    { 
	    *newcycle = TRUE;	
	    *p	= minp;
	    *c 	= minc; 
	    *o 	= mino;	

        if (minp <= 2)
        {
            matparmdone = FALSE;
        }
        else
        {
            matparmdone = TRUE;
        }

	    started = TRUE;	
	}
	else 
    {	
        /* ask user if we should continue */
	    getout=ui_estimate(*p,*c,*o,0); 
	    if (getout)	
        { 
		    *o = getobj;
		    return(TRUE);
	    }
   
        /* calculate the next set of indices */
        (*c)++;
        if (*c>maxc)
        {
            *c = minc;
            (*p)++;
        }
        if ((*p > 2) && (*o < nobject) && (*o <= maxo) && (!matparmdone))
        {
            *p = minp;
            (*o)++;
        }
        while ((!matparmdone) && ((!object[*o]) || (*o == ELMO)) && (*o < nobject) && (*o <= maxo)) 
               (*o)++;

        if ((*o == nobject) || (*o == maxo))
        {
             matparmdone = TRUE;
             *o = 0;
             *p = 3;
             *c = 0;
        }
        if (matparmdone)
        {
            if (*p > maxp) *newcycle = TRUE;
        }      

	    if (*newcycle)
        {
            cycle++;
	    	/*  check with user */
	    	getout=ui_estimate(*p, *c, *o, 1); 
	        if (getout)
            { /*reset */
		    	*o = getobj;
	    		return(TRUE);
	        }

		    /* restart */
	        *p	= minp;
	        *c 	= minc; 

            if (minp <= 2)
            {
	            *o 	= mino;	
                matparmdone = FALSE;
            }
            else
            {
                *o = maxo;
                matparmdone = TRUE;
            }


	    }
	}
	/* tell the parameter value */
	ui_estimate(*p, *c, *o, 2); 

	return(FALSE);  /* means:  don't stop */
}
/*------------------------------------------------------------------------*/
int changeestimate(
                   int ptype,
                   int pcoef,
                   int obj,
                   int st)
/*------------------------------------------------------------------------*/
{
int done;
int co1, co2;


/* perform the estimation and set the currentparm global */
/* when newparm is true, it tells the estimator to create a new stack */
currentparm = estimator(ptype, obj, pcoef, newparm, st, &done);

/* reset the actual parameters to the new parameter value */
getcoupledobjs(ptype,&co1,&co2);
switch(ptype)
{
    case MASS: 
    { 
    	mass[obj][pcoef] = GetPlateMaterial(ptype,pcoef,currentparm,obj); 
		break; 
    }  
    case RESIST: 
    { 
    	resist[obj][pcoef] = GetPlateMaterial(ptype,pcoef,currentparm,obj); 
		break; 
    }
    case STIFF: 
    { 
    	stiff[obj][pcoef] = GetPlateMaterial(ptype,pcoef,currentparm,obj); 
		break; 
    } 
    case cACRLTMstiff: 
    case cACRLOHCstiff: 
    case cACRLBAMstiff: 
    case cACRLACRLstiff:
    case cOHCOHCstiff: 
    case cOHCBAMstiff: 
    { 
        cstiff[co1][co2][pcoef]  = cstiff[co2][co1][pcoef]  = currentparm; 
        break;
    }
    case cACRLTMresist: 
    case cOHCBAMresist: 
    { 
        cresist[co1][co2][pcoef] = cresist[co2][co1][pcoef] = currentparm; 
        break;
    }
	default: {error("Parameter type is undefined.","in changeestimate()");}
}

return(done);

}
/*------------------------------------------------------------------------*/
double estimator(
                 int ptype, 
                 int obj,
                 int pcoef, 
                 int newstack, 
                 int stability, 
                 int *done) 
/*------------------------------------------------------------------------*/
{
    static double high,low;
    static int cla,clb;
    static int lowbounded, hibounded;
    static int started, i,mn, nstable;
    static struct stacktype stack[nestimates];
    double BIG_ERROR, newp, mindiff;
    double meshn[nestimates];
    double a,aa,b,bb,min;
    int totn;
    int rightbracket,leftbracket;
    int lb,rb,cp,j,cpfound;
    int equals;

    *done=FALSE;
    BIG_ERROR=1000.0;
    if (newstack) { i=0; mn=0; started=FALSE;}
    totn=meshsize = 2;


    /* initialize the boundaries */
    if (i==0) {
        if (ptype > 2) 
            obj = nobject;
	    high=range[hi][ptype][obj][pcoef];
	    low=range[lo][ptype][obj][pcoef];
	    if (low == high) *done = TRUE;
	    }

    for (j=0;j<=totn;j++) {
        meshn[j]=low+j*((high-low)/totn);
        if (!started) {
	        if (j==0) fprintf(stderr,"\n");
	        fprintf(stderr,"\t meshn[%d]=%e \n",j,meshn[j]); 
	        if (j==totn) fprintf(stderr,"\n");
	        fprintf(estimatefile.fp,"\t meshn[%d]=%e \n",j,meshn[j]); 
	        switch(ptype) 
            {
	            case MASS: 
                {
		            mass[obj][pcoef]=GetPlateMaterial(ptype,pcoef,meshn[j],obj); 
		            break;
	            }
	            case RESIST: 
                {
                    resist[obj][pcoef]=GetPlateMaterial(ptype,pcoef,meshn[j],obj); 
                    break;
	            }
	            case STIFF: 
                {
	 	            stiff[obj][pcoef]=GetPlateMaterial(ptype,pcoef,meshn[j],obj);
	 	            break;
	            }
	            default: { /* do nothing */ }
	        }	
        } 
    }
    mindiff=(high-low)/20;

    if (started) /* set up stack */
    { 
        stack[i].p = currentparm;

        if (nmodel>0) 
            stack[i].x = success[nscanmodel-1][nmodel-1].scan_norm_tot;
    
        if (stability != stable) 
        {
	        stack[i].x = BIG_ERROR;
	        nstable++;
	        /* find a likely stable estimate nearby */
	        if (mn <= totn) 
            {   /* if the mesh isn't explored, go to the next value */ 
	            if (stack[i].p != high) 
                    newp = meshn[mn]; 
	        }
	        else 
            {
	            min     = high - low;
	            cpfound = FALSE;
	            for (j = 0;j < i;j++) 
                {
	               if ((stack[j].x != BIG_ERROR) 
                    && (fabs(stack[i].p-stack[j].p) < min)) 
                   {
			           min      = fabs(stack[i].p - stack[j].p);
			           cp       = j;
			           cpfound  = TRUE;
			       }
		        }
		        if (cpfound) 
                {
			        newp=stack[i].p-(stack[i].p-stack[cp].p)/2.0;
		        }
		        else 
                {
		            /* give up, the source of the instability is another
		               parameter.  So put this parameter back to something
		               close to its initial value */
		            *done=TRUE;
		            newp = resetparm(ptype,pcoef,obj);
		            return(newp);
		        }
	        }
	        /* and get out of here  - and don't increment the stack */
	        mn++;
	        return(newp);
	    } /* end of instability handling */

        /* sort the stack so the the top of the stack is the lowest x */
        equals= sort(stack,i);
        /* make sure its not all the same */
        if ((equals == 1) && (i == totn)) 
        {
	        *done=TRUE;
	        newp=resetparm(ptype,pcoef,obj);
	        return(newp);
	    }
        /* increment the stack pointer */
        i++;
    }
    else started = TRUE;

    /* do the mesh parameters without any ado */
    if (mn<=totn) {
	    newp=meshn[mn];
	    mn++;
	    return(newp);
	    }

    /* first, find the p's closest to, and bracketing, the best current p */
    aa=bb=high-low;
    rightbracket=leftbracket=NOTFOUND;
    lb=rb=-1;
    for (j=0;j<i;j++) 
    {
	    if (stack[0].p<stack[j].p) 
        {
		    a=stack[j].p - stack[0].p;
		    rightbracket=FOUND;
		    if (a<aa) { aa = a; rb = j; }
	    }
	    if (stack[0].p>stack[j].p) 
        {
		    b=stack[0].p - stack[j].p;
		    leftbracket=FOUND;
		    if (b<bb) {bb = b; lb = j; }
	    }
    }	

    /* case 0 - both rb and lb are -1 so the last step was larger than 
    the margin=high-low... */
    if ((rb==-1) && (lb==-1)) 
    {
	    *done=TRUE;
	    newp=stack[0].p;
    }
    /* case 1  - there is no left bracket, p is smallest value in stack */
    else if (leftbracket==NOTFOUND)
    {
	    newp=stack[0].p-aa/2.0;
	    cla=clb=0;
    }
    /* case 2  - there is no right bracket, p is largest value in stack */
    else if (rightbracket==NOTFOUND) 
    {
	    newp=stack[0].p+bb/2.0;
	    cla=clb=0;
    }
    /* case 3  - p is bracketed by p[lb] and p[rb] */
    /* bracket on the right is closest - index of 1 - and only 2 attempts */
    else if ((rb<lb)&&(clb<2)) 
    {
	    /* put newp between p[rb] and p[0], closer to p[0] */ 
	    newp=stack[0].p + (stack[rb].p-stack[0].p)/3.0;
	    clb++;
	    cla=0;
    }
    /* bracket on the left is closest - index of 1 - and only 2 attempts */
    else if ((lb<rb)&&(cla<2)) 
    {
	    /* put newp between p[0] and p[lb], closer to p[0] */ 
	    newp=stack[0].p - (stack[0].p-stack[lb].p)/3.0;
	    cla++;
	    clb=0;
    }
    /* bracket on the left is closest but 2 attempts have been made on left */
    else if ((lb<rb)&&(cla>=2)) 
    { 
	    /* put newp between p[0] and p[rb], closer to p[0] */ 
	    newp=stack[0].p + (stack[rb].p-stack[0].p)/3.0;
	    cla=0;
    }
    /* bracket on the right is closest but 2 attempts have been made on right */
    else if ((rb<lb)&&(clb>=2)) 
    { 
	    /* put newp between p[0] and p[lb], closer to p[0] */ 
	    newp=stack[0].p - (stack[0].p-stack[lb].p)/3.0;
	    clb=0;
    }
    /* if best p is at a boundary, go back a bit and make sure its really minimal */
    if (stack[0].p>=high)  
    {
	    /* 0 is the right bracket, newp is a little less than p[0] */
	    if (!hibounded) 
        {
		    newp=stack[0].p - (stack[0].p-stack[lb].p)/6.0;
		    cla=clb=0;
		    hibounded=TRUE;
	    }
	    else 
        {
		    *done=TRUE;
		    newp=stack[0].p;
	    }
    }
    else if (stack[0].p<=low)  
    {
	    /* 0 is the left bracket, newp is a little larger than p[0] */
	    if (!lowbounded) 
        {
		    newp=stack[0].p + (stack[rb].p-stack[0].p)/6.0;
		    cla=clb=0;
		    lowbounded=TRUE;
        }
	    else 
        {
		    *done=TRUE;
		    newp=stack[0].p;	
        }
	    
    }
    else 
    {
	    hibounded=FALSE;
	    lowbounded=FALSE;
    }

    if (fabs(stack[1].p-stack[0].p)<mindiff) 
    {
	    *done=TRUE;
	    newp=stack[0].p;
    }

    return(newp);
}
/*------------------------------------------------------------------------*/
int sort(
         struct stacktype s[nestimates],
         int n)
/*------------------------------------------------------------------------*/
{
int i,j;
struct stacktype hold;
int count;
FILE *fp;

for (i=0;i<=n;i++) {
    for (j=i;j<=n;j++) {
	if (s[j].x<s[i].x) {
	    hold=s[i];
	    s[i]=s[j];
	    s[j]=hold;
	    }
	}
    }

count=0;
for (i=0;i<=n;i++) {
    if (s[0].x==s[i].x) {
	count++;
	if (count==n) return(1);
	}
    }

if (n>0) {
	fp=stderr;
	fprintf(fp,"\nStack \t Parameter \t Success ");
	for (j=0;j<i;j++) fprintf(fp,"\n %d \t %e \t %e ",j,s[j].p,s[j].x);
	fprintf(fp,"\n");
	fp=estimatefile.fp;
	fprintf(fp,"\nStack \t Parameter \t Success ");
	for (j=0;j<i;j++) fprintf(fp,"\n %d \t %e \t %e ",j,s[j].p,s[j].x);
	fprintf(fp,"\n");
	}

return(0);
}
/*-------------------------------------------------------------------------*/
/* the exps are
	M1 = H1 - beta1
	S1 = 3T1 - 5beta1
	R1 = 1/2(M1+S1) + d1

the amps are 
	M0= pi^2 b rho H0 / 8 beta0
	S0=pi^6 b E T0^3 / 96(1-nu^2) beta^5
	R0=sqrt(M0,S0)*d0

H : 	area OC/beta (thickness of plate representing cochlear partition)
beta: 	membrane width
T1:   	basilar membrane thickness
d:	    damping
b :	    channel width
rho :	material density
E :	    Young's modulus
nu :	Poisson's ratio 

amps (conjugate)
	b 	rho-H0 	T-E-nu 	d0
coefs (conjugate)
	H1 	T1    d1

assume beta known

this makes it

m0	rho-H0 
s0 	T-E-nu
R0  	d0

m1 	H1 - mass load due to OC and fluid (should be nonlinear with input)
s1 	T1 - thickness of the membrane contributes to plate stiffness
r1 	d1

the go-betweens are beta (membrane widths), b (channel width) and 

*/
/*-------------------------------------------------------------------------*/
double GetPlateMaterial(
                        int ptype,
                        int pcoef,
                        double parameter,
                        int obj)
/*-------------------------------------------------------------------------*/
{
    double pi6, beta0, beta1, beta5;
    double rho_H0, H1, T_E_NU, T1; 
    double d0, d1, insq;
    double b0;
    double outparm;

    if (!plate[obj]) return(parameter);

    /* init */
    pi6=M_PI*M_PI*M_PI*M_PI*M_PI*M_PI;

    /* membrane widths */
    beta0=width[obj][0];
    beta1=width[obj][1];
    beta5=width[obj][0]*width[obj][0]*width[obj][0]*width[obj][0]*width[obj][0];

    /* basal channel width */
    b0=(height[UC[obj]][0]+height[LC[obj]][0])/2.0;

    switch(ptype) {
	    case MASS: {
		    if (pcoef==0) { 
			    rho_H0=parameter; 
			    outparm = (M_PI*M_PI*b0*rho_H0)/(8.0*beta0); }	
		    else { 
			    H1=parameter;
			    outparm = H1-beta1; }
		    break; }
	    case STIFF: {
		    if (pcoef==0) {
			    /* this is really t^3 *E/(1-nu^2) */
			    T_E_NU=parameter;
			    outparm = (pi6*b0*T_E_NU)/(96.0*beta5); }
		    else {
			    T1=parameter;
			    outparm = (3*T1)-beta1; }
		    break; }
	    case RESIST: {
		    if (pcoef==0) {
			    d0=parameter;
			    insq=mass[BAM][0]*stiff[BAM][0];
			    outparm = sqrt(insq) *d0; }
		    else {
			    d1=parameter;
			    outparm =.5*(mass[BAM][1]+stiff[BAM][1])+d1; }
		    break; }
	    default:{ }
	    }

    return(outparm);    
}
/*-------------------------------------------------------------------------*/
double SetPlateMaterial(
                        int ptype,
                        int pcoef,
                        double inparm,
                        int obj)
/*-------------------------------------------------------------------------*/
{
    double pi6, beta0, beta1, beta5;
    double parameter;
    double rho_H0, H1, T_E_NU, T1; 
    double d0, d1, insq;
    double b0;

    if (!plate[obj]) return(inparm);

    /* init */
    pi6=M_PI*M_PI*M_PI*M_PI*M_PI*M_PI;

    /* membrane widths */
    beta0=width[obj][0];
    beta1=width[obj][1];
    beta5=width[obj][0]*width[obj][0]*width[obj][0]*width[obj][0]*width[obj][0];

    /* basal channel width */
    b0=(height[UC[obj]][0]+height[LC[obj]][0])/2.0;

    switch(ptype) {
	    case MASS: {
		    if (pcoef==0) { 
			    rho_H0= (8.0*beta0*inparm)/(M_PI*M_PI*b0);
			    parameter=rho_H0; 
			    }	
		    else { 
			    H1=inparm+beta1;  
			    parameter=H1;
			    }
		    break;
		    }
	    case STIFF: {
		    if (pcoef==0) {
			    /* this is really t^3 *E/(1-nu^2) */
			    T_E_NU=(inparm*96.0*beta5)/(pi6*b0);
			    parameter=T_E_NU;
			    }
		    else {
			    T1=(inparm+beta1)/3.0;
			    parameter=T1;
			    }
		    break;
		    }
	    case RESIST: {
		    if (pcoef==0) {
			    insq=mass[BAM][0]*stiff[BAM][0];
			    d0=inparm/ sqrt(insq) ;
			    parameter=d0;
			    }
		    else {
			    d1=inparm-.5*(mass[BAM][1]+stiff[BAM][1]);
			    /* notice, this is the only initial value that
			       can become zero 
			    if (d1==0) d1=inparm/2.0; */
			    parameter=d1;
			    }
		    break;
		    }
	    default:{
		    error("No such case in SetPlateMaterial()","in buildestimate.c");
		    }
	    }
    return(parameter);
}
/*------------------------------------------------------------------------*/
void setplateranges()
/*------------------------------------------------------------------------*/
{
int obj,i,j;
double initparm;

/* you only do this once, at the beginning of the whole estimation */
for (obj = 0; obj < nobject; obj++) {
    if ((object[obj])&&(plate[obj])) {
		for (i=0;i<3;i++) {
    			for (j=0;j<2;j++) {
					initparm=SetPlateMaterial(i,j,range[lo][i][obj][j],obj);
					range[lo][i][obj][j]=initparm;
					initparm=SetPlateMaterial(i,j,range[hi][i][BAM][j],obj);
					range[hi][i][obj][j]=initparm;
				}
			}
		}
	}
}

/*------------------------------------------------------------------------*/
void printMOREstuff(
                    FILE *fp,
                    int ptype,
                    int pcoef, 
                    int obj)
/*------------------------------------------------------------------------*/
{

	fprintf(fp,"The current parameter is: \t %s",ui_returnname(ptype));
    if (ptype <= 2) fprintf(fp,"[%d][%d]",obj,pcoef);

    fprintf(fp," = %e \n",ui_returnvalue(ptype,pcoef,obj));

}
/*------------------------------------------------------------------------*/
double ui_returnvalue(
               int p,
               int c, 
               int o)
/*------------------------------------------------------------------------*/
{
    int co1, co2;

    getcoupledobjs(p,&co1,&co2);

    switch(p) 
    {
	    case RESIST:    return(resist[o][c]); 
	    case STIFF:     return(stiff[o][c]);
	    case MASS:      return(mass[o][c]);
        /* coupled parameters */
        case cACRLTMstiff:  
        case cACRLOHCstiff: 
        case cACRLBAMstiff: 
        case cACRLACRLstiff:
        case cOHCOHCstiff:  
        case cOHCBAMstiff:   return(cstiff[co1][co2][c]); 

        case cACRLTMresist: 
        case cOHCBAMresist:  return(cresist[co1][co2][c]);
	    default: 
        { 
            error("Parameter type is undefined.", "in ui_returnvalue()");
            return(0);
        }
    }
}
/*------------------------------------------------------------------------*/
char *ui_returnname(
              int ptype)
/*-------------------------------------------------------------------------*/
{


    switch(ptype) 
    {
	    case RESIST:            return("resist"); 
	    case STIFF:             return("stiff"); 
	    case MASS:              return("mass");
        case cACRLTMstiff:      return("ACRL_TMstiff");
        case cACRLOHCstiff:     return("ACRL_OHCstiff");
        case cACRLBAMstiff:     return("ACRL_BAMstiff");
        case cACRLACRLstiff:    return("ACRL_ACRLstiff");
        case cOHCOHCstiff:      return("OHC_OHCstiff");
        case cOHCBAMstiff:      return("OHC_BAMstiff"); 

        case cACRLTMresist:     return("ACRL_TMresist");
        case cOHCBAMresist:     return("OHC_BAMresist");
	    default: 
        { 
		    error("Parameter type is undefined.", "in nametype()"); 
            return("");
        }
	}
}

/*------------------------------------------------------------------------*/
int ui_estimate(
                int p,
                int c,
                int o,
                int ui_type)
/*------------------------------------------------------------------------*/
{
    int co1, co2;

    switch(ui_type) 
    {
        case 0: 
        {
    	    shortsummary(estimatefile.fp,o,p);
	        fprintf(stderr,"\nThe final parameters for this pass are: \n");
    	    shortsummary(stderr,o,p);
    	    fprintf(stderr,"\nCHANGING PARAMETERS\n"); 
    	    fprintf(estimatefile.fp,"\nCHANGING PARAMETERS\n"); 
    	    if (pauses==1) 
            {
	            if ((query("Continue (0) or Stop (1)?",0,1))==1)  return(TRUE);
	        }
	        break;
	    }
        case 1: 
        {
	        fprintf(stderr,"The current cycle of estimation is complete.\n");
	        fprintf(stderr,"\nThe next CYCLE is %d\n",cycle+1); 
	        fprintf(estimatefile.fp,"\nThe next CYCLE is %d\n",cycle+1); 
	        if (pauses==2) 
            {
	            if ((query("Continue (0) or Stop (1)?",0,1))==1)  return(TRUE);
	        }
	        break;
	    }
        case 2: 
        {
            if (p>2)
                getcoupledobjs(p,&co1,&co2);
            else co1 = o;
	        fprintf(stderr, "For object ");
	        nameit(co1,stderr);
            if (p>2)
            {
                fprintf(stderr," coupled to ");
                nameit(co2,stderr);
            }
	        fprintf(stderr,"the parameter is %d, coef is %d \n",p,c);

	        fprintf(estimatefile.fp, "For object ");
	        nameit(co1,estimatefile.fp);
            if (p>2)
            {
                fprintf(estimatefile.fp," coupled to ");
                nameit(co2,estimatefile.fp);
            }
	        fprintf(estimatefile.fp,"the parameter is %d, coef is %d \n",p,c);


	    }
    } 
    return(FALSE);	
}
/*------------------------------------------------------------------------*/
void getbounds(
        int *mino,
        int *maxo,
        int *minp,
        int *maxp,
        int *minc,
        int *maxc)
/*------------------------------------------------------------------------*/
{
    *mino = estimator_bound[0];
    *maxo = estimator_bound[1];
    *minp = estimator_bound[2];
    *maxp = estimator_bound[3];
    *minc = estimator_bound[4];
    *maxc = estimator_bound[5];

}
/*------------------------------------------------------------------------*/
void setbounds(
        int mino,
        int maxo,
        int minp,
        int maxp,
        int minc,
        int maxc)
/*------------------------------------------------------------------------*/
{

    estimator_bound[0] = mino;
    estimator_bound[1] = maxo;
    estimator_bound[2] = minp;
    estimator_bound[3] = maxp;
    estimator_bound[4] = minc;
    estimator_bound[5] = maxc;

}