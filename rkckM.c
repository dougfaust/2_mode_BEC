
#include "math.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

#define NRANSI
#include "nrutil.h"

void rkckM(double y[], double dydx[], int n, double x, double h, double yout[],double yerr[], 
	void (*derivs2)(double [], double [], int) )
{	
	int i;	 
	static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
		b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
		b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
		b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
		b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
		c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
		dc5 = -277.00/14336.0;
	double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
		dc4=c4-13525.0/55296.0,dc6=c6-0.25;
	double *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;
	
	double re_err,im_err;

	ak2=vector(0,(long)n);
	ak3=vector(0,(long)n);
	ak4=vector(0,(long)n);
	ak5=vector(0,(long)n);
	ak6=vector(0,(long)n);
	
	ytemp=vector(0,(long)n);
	
	for (i=0;i<=n;i++)
		{	ytemp[i]=y[i]+b21*h*dydx[i];}
	
	(*derivs2)(ytemp, ak2, n);
	for (i=0;i<=n;i++)
		{ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);}
	
	(*derivs2)(ytemp, ak3, n);
	for (i=0;i<=n;i++)
		{ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);}
		
	(*derivs2)(ytemp, ak4, n);
	for (i=0;i<=n;i++)
		{ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);}
	
	(*derivs2)(ytemp, ak5, n);
	for (i=0;i<=n;i++)
		{ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);}
	
	(*derivs2)(ytemp, ak6, n);
	for (i=0;i<=n;i++)
		{yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);}
		
	for (i=0;i<=(n-1)/2;i++)
		{re_err=h*(dc1*dydx[2*(i-1)+1]+dc3*ak3[2*(i-1)+1]+dc4*ak4[2*(i-1)+1]+dc5*ak5[2*(i-1)+1]+dc6*ak6[2*(i-1)+1]);
		 im_err=h*(dc1*dydx[2*i]+dc3*ak3[2*i]+dc4*ak4[2*i]+dc5*ak5[2*i]+dc6*ak6[2*i]);
		yerr[i]=sqrt(re_err*re_err + im_err*im_err);
		}
	free_vector(ytemp,0,(long)n);
	free_vector(ak6,0,(long)n);
	free_vector(ak5,0,(long)n);
	free_vector(ak4,0,(long)n);
	free_vector(ak3,0,(long)n);
	free_vector(ak2,0,(long)n);
	
}
#undef NRANSI
