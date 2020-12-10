/* CAUTION: This is the ANSI C (only) version of the Numerical Recipes
   utility file nrutil.h.  Do not confuse this file with the same-named
   file nrutil.h that is supplied in the 'misc' subdirectory.
   *That* file is the one from the book, and contains both ANSI and
   traditional K&R versions, along with #ifdef macros to select the
   correct version.  *This* file contains only ANSI C.               */

#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))

static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static double minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

static long lmaxarg1,lmaxarg2;
#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
        (lmaxarg1) : (lmaxarg2))

static long lminarg1,lminarg2;
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
        (lminarg1) : (lminarg2))

static long imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

static long iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))




// complex (ie. actual) fft header

void four1(double data[], unsigned long nn, int isign);
void fourn(double data[], unsigned long nn[], int ndim,  int isign);
void sinft(double y[], int n);
void realft(double data[], unsigned long n, int isign);
void cfft_d_(double*,int*,int*);  //hand coded fast fourier transforms


// add headers for specific routines added from Num Rec nr.h

// wavefunction drivers

void odeintW(double ystart[], int nvarr, double x1, double x2, double eps, double h1, double hmin, int *nok, int *nbad, double ke_x[], 
                double scale_factor, unsigned long nvar_x, void (*derivs)(double,double [], double [], double [], double, unsigned long), 
                void (*rkqsW)(double [], double [], int, double *, double, double, double [], double *, double *, double[],double,unsigned long, 
                void (*)(double,double [],double [],double[],double,unsigned long)));

void rkckW(double y[], double dydx[], int n,double x, double h, double yout[],
	double yerr[], double ke_x[], double scale_factor, unsigned long nvar_x,
	void (*derivs)( double, double [], double [],double[],double,unsigned long) );
	
void rkqsW(double y[],double dydx[],int n,double *x,double htry
		  ,double eps, double yscal[],double *hdid, double *hnext,
		  double ke_x[], double scale_factor, unsigned long nvar_x,
		  void(*derivs)(double, double [],double [],double[],double,unsigned long));

void derivs(double x,double y[],double dydx[], double ke_x[], 
			double scale_factor, unsigned long nvar_x);

// CI drivers

void odeintM(double [], int, double, double, double, double, double, int *, int *, void (*derivs2)(double [], double [], int), 
		     void (*rkqsM)(double [], double [], int, double *, double, double, double [], double *, double *,
			  void (*)(double [], double [], int)));

void rkqsM(double y[],double dydx[],int n,double *x,double htry, double eps, double yscal[], double *hdid, double *hnext, void(*derivs2)(double [], double [], int));

void rkckM(double y[], double dydx[], int n, double x, double h, double yout[], double yerr[], void (*derivs2)(double [], double [], int));

void derivs2(double c[], double dcdt[], int n);	

// stiff equation routine headers.

void simpr(double y[],double dydx[],double dfdx[],double **dfdy, long n, double xs, double htot, 
		   long step, double yout[],void((*derivs)(double,double[],double[]))); 
void stifbs(double y[], double dydx[], int nv, double  *xx,
			double htry, double eps,double yscale[],double *hdid, double *hnext,
			void((*derivs)(double,double[],double[]))); 
void pzextr(long iest, double xest,double yest[],double yz[],double dy[],int nv);
void ludcmp( double **a,long n,long *indx,double *d);
void lubskp( double **a,long n,long *indx,double b[]);
void jacobn( double x, double y[], double dfdx[], double ** dfdy, long n);

double sech(double q);
void nrerror(char error_text[]);
double *vector(long nl, long nh);
int *ivector(long nl, long nh);
unsigned char *cvector(long nl, long nh);
unsigned long *lvector(long nl, long nh);
double *dvector(long nl, long nh);
double **matrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
long **imatrix(long nrl, long nrh, long ncl, long nch);
double **submatrix(double **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl);
double **convert_matrix(double *a, long nrl, long nrh, long ncl, long nch);
double ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_vector(double *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_cvector(unsigned char *v, long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(long **m, long nrl, long nrh, long ncl, long nch);
void free_submatrix(double **b, long nrl, long nrh, long ncl, long nch);
void free_convert_matrix(double **b, long nrl, long nrh, long ncl, long nch);
void free_f3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh);


#endif /* _NR_UTILS_H_ */
