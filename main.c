
/* 
   computational model of Bose-Einstein splitting experiments (c.f. https://doi.org/10.1103/PhysRevLett.99.240406)
   
   theory behind code detailed in https://arxiv.org/abs/1008.0217
   
   this version includes:
   rewrite of full WF and CI 2 mode dynamics with RK4 drivers for either variable, 
   sensible units (AMU, msec, microns), sensible variable names, sensible organization 
   
   exists as DKF's skeleton code for CS to port to full 2D and 3D
*/

#include "math.h"
#include <stdio.h>
#include "nrutil.h"
#include <stdlib.h>

#include <time.h>

double N, g, g12, L, C11, C12, E0, scale_factor;
double **H, *ke;

int main()
{

  // start timer
  clock_t ticks1, ticks2;
  ticks1 = clock();

	void odeintW(double ystart[], int nvarr, double x1, double x2, double eps, double h1, double hmin, int *nok, int *nbad, double ke_x[], 
                double scale_factor, unsigned long nvar_x, void (*derivs)(double,double [], double [], double [], double, unsigned long), 
                void (*rkqsW)(double [], double [], int, double *, double, double, double [], double *, double *, double[],double,unsigned long, 
                void (*)(double,double [],double [],double[],double,unsigned long)));

	void sinft(double [], int);	
	void cxxcoeffs(double [], int);

        void odeintM(double [], int, double, double, double, double, double, int *, int *, void (*derivs2)(double [], double [], int), 
		     void (*rkqsM)(double [], double [], int, double *, double, double, double [], double *, double *,
			  void (*)(double [], double [], int)));

	void fillHam(double [], double, int, int);


	int i, nsteps, n, nok, nbad, np, nv;
	double h1, hmin, timestep, eps, temp;
	double *c, *a, *a_i, *a_r;
	double norm, over_R, over_I, time, dx, xx, norm_factor;
	double Pi = 3.141592653589793;
	double Csec, Ssec, temp1, temp2;

	char *wR = "WFRealPart.out";
	char *wI = "WFImgPart.out";
	char *wD = "WF_Density.out";
	char *inW = "InputWF.dat";
	char *cR = "CIRealPart.out";
	char *cI = "CIImgPart.out";
	char *cD = "CI_Density.out";
	char *inC = "InputCI.dat";
	char line[22];

	FILE *fcR, *fcI, *fcD, *readC;
	FILE *fwR, *fwI, *fwD, *readW;

	extern double scale_factor, g, g12, L, N, C12, C11, E0, **H, *ke;
	unsigned long npart, nvar_x;

	nvar_x = 256;                        // grid resolution
	nv = (int)nvar_x;
	L = 25.0;                            // boxlength, in microns
	dx = L/(double)nv;
	norm_factor =  (2.0/(double)(nv));
	scale_factor = 2.76253;              // this is h_bar/M for 23Na

	N = 2.0e2;                           // number of particles
	npart = (unsigned long)N;
	np = (int)N;
	fprintf(stdout, "\nN = %f,  npart = %u, np = %d\n", N, npart, np);
	     
	a = vector((long)1, 2*nvar_x);
	a_i = vector((long)1, nvar_x);
	a_r = vector((long)1, nvar_x);
	ke = vector((long)1, nvar_x);
	c = vector((long)0, 2*npart+1);  // fock space vector, CI
	H = matrix((long)0, 2*npart+1, (long)0, 2*npart+1);  // fock space Hamiltonian

	nsteps = 550;
	timestep = 0.2;
	h1 = 0.001;
	hmin = 1.0e-12;
	eps = 2.0e-5;
	g = 8.06496e-3;  // this is usual 'g' over the transverse TF radii 
    g12 = 1.2*g;

	E0 = 22500.0;  // this should be manually entered.  H[N/2][N/2] = E0 is a good starting guess
	// E0 = 0.000;

	Csec = cos(E0*timestep/63.5077);
	Ssec = sin(-E0*timestep/63.5077);
	fprintf(stdout, " energy offset = %f \n", E0);

	for(i=1;i<=nv;i++)  ke[i] = scale_factor*((1.0-cos((double)(i-1)*Pi/(double)(nv)))/(dx*dx));


	// Open file for write and output initial vars
	fprintf ( stdout, "opening files for read/write \n");
	readC=fopen(inC, "r");	
	fcD=fopen(cD, "w");
	fcI=fopen(cI, "w");
	fcR=fopen(cR, "w");
	readW=fopen(inW, "r");	
	fwD=fopen(wD, "w");
	fwI=fopen(wI, "w");
	fwR=fopen(wR, "w");

	// read in initial WF from file
        for(i=2;i<=2*nv;i+=2)
        {
	 fgets(line, 22, readW);
	 a[i-1]=atof(line);
	 a[i]=0.00;
	 a_r[i/2] = a[i-1];
	 a_i[i/2] = a[i];
        }	

	// read in initial CI distribution from file
        for(i=0;i<=(2*np+1);i+=2)
        {
	 fgets(line, 22, readC);
	 c[i]=atof(line);
	 c[i+1]=0.00;
        }
	
	// check that initial CI dist is a properly-normalized thing
	norm = 0.00;
	for(i=0;i<=2*np+1;i++) norm += c[i]*c[i];
	norm = sqrt(norm);
	fprintf(stdout, "CI norm = %e   ", norm);

	// check that initial WF is a properly-normalized thing, orthogonal to its reflection o'er the center of the trap
	norm = 0.00;
	over_R = 0.00;
	over_I = 0.00;
	for(i=1;i<=2*nv;i+=2)
	  {
	    norm += (a[i]*a[i] + a[i+1]*a[i+1]);
	    over_R += (a[2*nv-i]*a[i] + a[2*nv-i+1]*a[i+1]);
	    over_I += (a[2*nv-i]*a[i+1] - a[2*nv-i+1]*a[i]);
	  }
	over_R *= dx; over_I *= dx; norm *= dx;
	norm = sqrt(norm);
	fprintf(stdout, "WF norm = %e\n WF overlap = (%e) + i(%e)\n", norm, over_R, over_I);

	// debuggery
	for(i=1;i<=2*nv;i++) a[i] /= norm;
	norm = 0.00;
	over_R = 0.00;
	over_I = 0.00;
	for(i=1;i<=2*nv;i+=2)
	  {
	    norm += (a[i]*a[i] + a[i+1]*a[i+1]);
	    over_R += (a[2*nv-i]*a[i] + a[2*nv-i+1]*a[i+1]);
	    over_I += (a[2*nv-i]*a[i+1] - a[2*nv-i+1]*a[i]);
	  }
	norm = sqrt(norm);
	fprintf(stdout, "\n\n NEW WF norm = %e\n WF overlap = (%e) + i(%e)\n", norm, over_R, over_I);

	// compute CXX coefficients
	cxxcoeffs(c, np);
	fprintf(stdout, "C11 = %20.15f  C12 = %20.15f \n", C11, C12);

	/* DON'T print out initial state 
	xx = 0.00;
        for(i=1;i<=2*nv;i+=2)
        {
	  fprintf (fwR, "0.0 %f  %20.8f\n", xx, a[i]);  // real part of WF
	  fprintf (fwI, "0.0 %f  %20.8f\n", xx, a[i+1]);  // img part of WF
	  fprintf(fwD, "0.0 %f %20.8f\n", xx, (C11*(a[i+1]*a[i+1]+a[i]*a[i]) + C12*(a[i]*a[2*nv-i]+a[i+1]*a[2*nv-i+1]) + (N-C11)*(a[2*nv-i+1]*a[2*nv-i+1]+a[2*nv-i]*a[2*nv-i]))); // density
	  xx += dx;
        }	 
*/

	// put WF into k-space for faster convergence
	sinft(a_r, nv);
	sinft(a_i, nv);
	for(i=2;i<=2*nv;i+=2)
	{
		a[i]   = a_i[i/2];
		a[i-1] = a_r[i/2];
	}


	/* ODEINT LOOP TO PROPAGATE STATE TO FINAL TIME */
	time = 0.00;
	for(n=1;n<=nsteps;n++)
	  {

	    // take a step for the WF and return to x-space
	    (*odeintW)(a,  2*nv, time, (time+timestep), eps,  h1, hmin,  &nok,  &nbad, ke, scale_factor, nvar_x, derivs,  rkqsW);
		for(i=2;i<=2*nv;i+=2)
		{
		  a_i[i/2] = norm_factor*a[i];
		  a_r[i/2] = norm_factor*a[i-1];
		}
		sinft(a_r, nv);
		sinft(a_i, nv);
		for(i=2;i<=2*nv;i+=2)
		{
		  a[i] = a_i[i/2];
		  a[i-1] = a_r[i/2];
		}
	        // check WF norm and overlap
		norm = 0.00;
		over_R = 0.00;
		over_I = 0.00;
		for(i=1;i<=2*nv;i+=2)
		  {
		    norm += (a[i]*a[i] + a[i+1]*a[i+1]);
		    over_R += (a[2*nv-i]*a[i] + a[2*nv-i+1]*a[i+1]);
		    over_I += (a[2*nv-i]*a[i+1] - a[2*nv-i+1]*a[i]);
		  }
		over_R *= dx; over_I *= dx; norm *= dx;
		norm = sqrt(norm);
		fprintf(stdout, "WF norm = %e\n WF overlap = (%e) + i(%e)\n", norm, over_R, over_I);


	    // take a step for the CI
	    fillHam(a, time, np, nv);
	    (*odeintM)(c, (2*np+1), time, (time+timestep), eps, h1, hmin, &nok, &nbad, derivs2, rkqsM);
	           // check CI norm
	           norm = 0.00;
	           for(i=0;i<=2*npart+1;i++) norm += c[i]*c[i];
	           norm = sqrt(norm);
	           fprintf(stdout, "CI norm = %e\n", norm);
               for(i=0;i<=2*np+1;i++) c[i] /= norm;


	    /* put in secular phase evolution for large number of particles this is NOT recommended
	    for(i=0;i<=2*np+1;i+=2)
	      {
		temp1 = Csec*c[i]-Ssec*c[i+1];
		temp2 = Ssec*c[i]+Csec*c[i+1];
		c[i] = temp1;
		c[i+1] = temp2;
	      } */

	    time += timestep;
	    cxxcoeffs(c, np);

	    fprintf(stdout, "time = %f millisec, printing to file \n", time);

//
// skip the first 10 ms (added 5/27/2010 to show PRL referees that we do indeed have the true ground state)
//   
if(time >= 10.00){
	    // print output to appropriate files
	    fprintf(fcD, " \n ");
	    fprintf(fcI, " \n ");
	    fprintf(fcR, " \n ");
	    for(i=0;i<=2*np+1;i+=2)
	      { 
		fprintf(fcD, "%f  %f  %20.12f\n", (time-10.0), (float)i/2, c[i]*c[i]+c[i+1]*c[i+1]);
		fprintf(fcI, "%f  %f  %20.12f\n", (time-10.0), (float)i/2, c[i+1]);
		fprintf(fcR, "%f  %f  %20.12f\n", (time-10.0), (float)i/2, c[i]);
	      }
	    fprintf(fwD, "\n");
	    fprintf(fwI, "\n");
	    fprintf(fwR, "\n");	    
	    xx = 0.00;
	    for(i=1;i<=2*nv;i+=2)
	      {
		fprintf(fwR, "%f  %f  %20.8f\n", (time-10.0), xx, a[i]);
		fprintf(fwI, "%f  %f  %20.8f\n", (time-10.0), xx, a[i+1]);
		fprintf(fwD, "%f  %f  %20.8f\n", (time-10.0), xx, fabs((C11*(a[i+1]*a[i+1]+a[i]*a[i]) + C12*(a[i]*a[2*nv-i]+a[i+1]*a[2*nv-i+1]) + (N-C11)*(a[2*nv-i+1]*a[2*nv-i+1]+a[2*nv-i]*a[2*nv-i]))));
		xx += dx;
	      }
}

	    // return WF to k-space for continued propagation if it's not the last timestep
		
	    if(n == nsteps) break;

	    for(i=2;i<=2*nv;i+=2)
	      {
		a_i[i/2] = a[i];
		a_r[i/2] = a[i-1];
	      }
            // return to k-space
	    sinft(a_r, nv);
	    sinft(a_i, nv);
	    for(i=2;i<=2*nv;i+=2)
	      {
		a[i]   = a_i[i/2];
		a[i-1] = a_r[i/2];
	      }



	  }


	fclose(fwD);
	fclose(fwR);
	fclose(fwI);
	fclose(readW);
	fclose(fcD);
	fclose(fcR);
	fclose(fcI);
	fclose(readC);
	fprintf(stdout, "file closed!\n");
	free_vector(a, 1, 2*nvar_x);
	free_vector(a_r, 1, nvar_x);
	free_vector(a_i, 1, nvar_x);
	free_vector(ke, 1, nvar_x);
	free_vector(c, 0, 2*npart+1);
	free_matrix(H, 0, 2*npart+1, 0, 2*npart+1);

	// stop clock
	ticks2 = clock();
	printf("\n\nTook %ld ticks to complete the program.\n",ticks2-ticks1);
	printf("this corresponds to %ld seconds.\n\n", (ticks2-ticks1)/CLOCKS_PER_SEC);

	return 0;

}

void cxxcoeffs(double c[], int np)  // these are the coefficients needed for the physical density
{
  extern double C11, C12;
  int j;
  C11 = 0.00;
  for(j=0;j<=2*np+1;j+=2){ C11 += (j/2.0)*(c[j]*c[j]+c[j+1]*c[j+1]); }
  C12 = 0.00;
  for(j=0;j<=2*np+1;j+=2){ C12 += sqrt(j/2.0+1.0)*sqrt(np-j/2.0)*(c[j]*c[j+2]+c[j+1]*c[j+3]); }
  C12 *= 2.0;
}

void derivs2(double c[], double dcdt[], int n)  // this is the Fock-space Hamiltonian part dcdt = -iH*c/hbar
{

  extern double **H;
  double temp, hbar = 63.5077; // AMU, millisec, microns
  int i, j;

  for(i=0;i<=n;i++){
    temp = 0.00;
    for(j=0;j<=n;j++)
      {
	temp += H[i][j]*c[j];
      }
    if(i%2 == 0){ dcdt[i+1] = -1.0*temp/hbar; } // fprintf(stdout, "filled dcdt[%d] = %e \n", i+1, dcdt[i+1]);}
    if(i%2 == 1){ dcdt[i-1] = temp/hbar; } // fprintf(stdout, "filled dcdt[%d] = %e \n", i-1, dcdt[i-1]);}
    }

  /*
  fprintf(stdout, "  \n ");
    fprintf(stdout, " derivs gives dcdt[i] = \n");
  for(i=0;i<=n;i++) fprintf(stdout, "  %e  ", dcdt[i]);
  */

}

void derivs(double t,double y[],double dydx[], double ke_x[], double scale_factor, unsigned long nvar_x)   // derivs for the wavefunction component
{
        int i,xind,index, nv; 
	double *a_r,*a_i,b,xx,mu,norm_factor,dx, norm;
	double *aTr, *aTi, rMu12, iMu12;

	extern double g, g12, L, N, *ke;

	void sinft(double [],int);
	double extpot(double, double);

	nv = (int)nvar_x;
	dx = L/((double)nvar_x);
	a_r = vector((long)1,(long)(nvar_x));
	a_i = vector((long)1,(long)(nvar_x));
	aTr = vector((long)1,(long)(nvar_x));
	aTi = vector((long)1,(long)(nvar_x));

	// enter in k-space
	norm_factor =  (2.0/(double)(nvar_x));
	for(i=2; i<=(2*nv); i+=2)
	{
	  a_i[i/2] = y[i]*norm_factor;
	  a_r[i/2] = y[i-1]*norm_factor;
	  aTi[i/2] = ke[i/2]*y[i]*norm_factor;
	  aTr[i/2] = ke[i/2]*y[i-1]*norm_factor;
	}
		
	sinft(a_r, nv);
	sinft(a_i, nv);
	sinft(aTr, nv);
	sinft(aTi, nv);

    norm = 0.00;
    for(i=1;i<=nv;i++) { norm += (a_i[i]*a_i[i]+a_r[i]*a_r[i]); }
    norm *= dx;
	//norm = sqrt(norm);

	xx = 0.00;
	for(i=1; i<=nv; i++)
	{
	  mu = (extpot(t, xx)/scale_factor);  
	  b = (g*(N/2.0-1.0)*(a_r[i]*a_r[i] + a_i[i]*a_i[i]) + g12*N*(a_r[nv-i+1]*a_r[nv-i+1] + a_i[nv-i+1]*a_i[nv-i+1]))*scale_factor/norm;  
	  aTr[i]  += (b + mu)*a_r[i];  
	  aTi[i] += (b + mu)*a_i[i];

	  //	  dydx[2*i-1] = a_r[i];  // deprecated way of storing derivs inherited from WPR's post-doc. i don't see why it's useful? 
	  //dydx[2*i] = a_i[i];		

	  xx += dx;
	}

        // Mu12 lagrange multiplier constraining the two orbitals to be orthogonal
	rMu12 = 0.00;
	iMu12 = 0.00;
	for(i=1;i<=nv;i++)
	{
	  rMu12 += (a_r[nv-i+1]*aTr[i] + a_i[nv-i+1]*aTi[i]);
	  iMu12 += (a_r[nv-i+1]*aTi[i] - a_i[nv-i+1]*aTr[i]);
	}
	rMu12 *= dx/norm;
	iMu12 *= dx/norm;

	for(i=1;i<=nv;i++)  
	{
	  aTr[i] += (iMu12*a_i[nv-i+1] - rMu12*a_r[nv-i+1]);
	  aTi[i] -= (rMu12*a_i[nv-i+1] + iMu12*a_r[nv-i+1]);
	}

	sinft(aTr, nv);
	sinft(aTi, nv);
	
	for(i=1;i<=nv;i++)
	{
	  dydx[2*i] = -1.0*aTr[i];
	  dydx[2*i-1] = aTi[i];
        }

	free_vector(a_r,(long)1,(long)(nvar_x));
	free_vector(a_i,(long)1,(long)(nvar_x));	
	free_vector(aTr,(long)1,(long)(nvar_x));
	free_vector(aTi,(long)1,(long)(nvar_x));
}

double extpot(double t, double x)
{

   extern double L;
   double KK, AMP, AA;

   KK = 13.9317;
   AMP = 1200.00; //1000.0/2.0*(1.0 + tanh(0.5*(t - 6.0)));
   AA = 0.1;
   //AMP = 0.00;

   return 0.5*KK*(x - 0.5*L)*(x - 0.5*L) + (1.0 + 0.1*(AMP))*exp(-AA*(x - 0.5*L)*(x - 0.5*L)) - 0.5*KK/AA*(1.0 - log(KK) + log(2.0*AA*(1.0 + 0.1*(AMP))));   
   // return 0.5*KK*(x-L/2.0)*(x-L/2.0);
   // return 0.00;  
}


void fillHam(double a[], double t, int np, int nvar)  // fill the CI basis Hamiltonian with wavefunction data
{

  extern double g, L, **H, E0;
  double G1111, G1221, rG1112, iG1112, rG1122, iG1122;  // T0, T1, T2 analogues, not real since we're TD
   double E11, rE12, iE12;
   double mu, dx, engys[3];
   int i;

   dx = L/(double)nvar;
   void fillEngys(double [], double [], double, int);


         // compute Gammas.  
	 // this code doesn't really look "D.R.Y.", but there are lots of non-trivial permutations between each GXXX.

  G1111 = 0.000;
 for(i=1;i<=2*nvar;i+=2)
   { 
     G1111 += (a[i]*a[i]*a[i]*a[i] + 2.0*a[i]*a[i]*a[i+1]*a[i+1] + a[i+1]*a[i+1]*a[i+1]*a[i+1]);
   }
 G1221 = 0.00;
 for(i=1;i<=2*nvar;i+=2)
   {
     G1221 += (a[i]*a[2*nvar-i]*a[2*nvar-i]*a[i] - a[i]*a[2*nvar-i]*a[2*nvar+1-i]*a[i+1] + a[i]*a[2*nvar+1-i]*a[2*nvar+1-i]*a[i] +
	       a[i]*a[2*nvar-i+1]*a[2*nvar-i]*a[i+1] + a[i+1]*a[2*nvar-i]*a[2*nvar-i+1]*a[i] + a[i+1]*a[2*nvar-i]*a[2*nvar-i]*a[i+1]- 
	       a[i+1]*a[2*nvar-i+1]*a[2*nvar-i]*a[i] + a[i+1]*a[2*nvar-i+1]*a[2*nvar-i+1]*a[i+1]);
   }
 rG1112 = 0.000;
 for(i=1;i<=2*nvar;i+=2)
   {
     rG1112 += (a[i]*a[i]*a[i]*a[2*nvar-i] - a[i]*a[i]*a[2*nvar+1-i]*a[i+1] + a[i]*a[i+1]*a[i+1]*a[2*nvar-i] +
	       a[i]*a[i+1]*a[i]*a[2*nvar-i+1] + a[i+1]*a[i]*a[i+1]*a[2*nvar-i] + a[i+1]*a[i]*a[i]*a[2*nvar-i+1] - 
	       a[i+1]*a[i+1]*a[i]*a[2*nvar-i] + a[i+1]*a[i+1]*a[2*nvar-i+1]*a[i+1]);
   }
 rG1122 = 0.000;
 for(i=1;i<=2*nvar;i+=2)
   {
     rG1122 += (a[i]*a[i]*a[2*nvar-i]*a[2*nvar-i] - a[i]*a[i]*a[2*nvar+1-i]*a[2*nvar-i+1] + a[i]*a[i+1]*a[2*nvar-i+1]*a[2*nvar-i] +
	       a[i]*a[i+1]*a[2*nvar-i]*a[2*nvar-i+1] + a[i+1]*a[i]*a[2*nvar-i+1]*a[2*nvar-i] + a[i+1]*a[i]*a[2*nvar-i]*a[2*nvar-i+1] - 
	       a[i+1]*a[i+1]*a[2*nvar-i]*a[2*nvar-i] + a[i+1]*a[i+1]*a[2*nvar-i+1]*a[2*nvar-i+1]);
   }

 iG1112 = 0.000;
 for(i=1;i<=2*nvar;i+=2)
   {
     iG1112 += (-a[i]*a[i]*a[i]*a[2*nvar-i+1] - a[i]*a[i]*a[2*nvar+1-i]*a[i] + a[i]*a[i+1]*a[i]*a[2*nvar-i] +
	       a[i+1]*a[i]*a[i]*a[2*nvar-i] + a[i+1]*a[i+1]*a[i+1]*a[2*nvar-i] + a[i+1]*a[i+1]*a[i]*a[2*nvar-i+1] - 
	       a[i+1]*a[i+1]*a[i]*a[2*nvar-i+1] - a[i]*a[i+1]*a[2*nvar-i+1]*a[i+1]);
   }
 iG1122 = 0.000;
 for(i=1;i<=2*nvar;i+=2)
   {
     iG1122 += (-a[i]*a[i]*a[2*nvar-i]*a[2*nvar-i+1] - a[i]*a[i]*a[2*nvar+1-i]*a[2*nvar-i] + a[i]*a[i+1]*a[2*nvar-i]*a[2*nvar-i] +
	       a[i+1]*a[i]*a[2*nvar-i]*a[2*nvar-i] + a[i+1]*a[i+1]*a[2*nvar-i+1]*a[2*nvar-i] + a[i+1]*a[i+1]*a[2*nvar-i]*a[2*nvar-i+1] - 
	       a[i+1]*a[i]*a[2*nvar-i+1]*a[2*nvar-i+1] - a[i]*a[i+1]*a[2*nvar-i+1]*a[2*nvar-i+1]);
   }

 G1111 *= dx; G1221 *= dx; rG1112 *= dx; iG1112 *= dx; rG1122 *= dx; iG1122 *= dx;
 
  // compute engys
 fillEngys(engys, a, t, nvar);
 E11 = engys[0];  
 rE12 = engys[1];
 iE12 = engys[2];
 // iE12 = 0.5; // DEBUGGERY
// fprintf(stdout, "\nthis is E11 %e\n", E11);
// fprintf(stdout, "this is E12 (%e)+i(%e)\n", rE12, iE12);


/* PHASE SHIFT */
double temp1, temp2;

if(t >= 0.00)
{
temp1 = cos(1.0*3.1415)*rE12-sin(1.0*3.1415)*iE12;
temp2 = sin(1.0*3.1415)*iE12+cos(1.0*3.1415)*rE12;
rE12 = temp1;
iE12 = temp2;

temp1 = cos(1.0*3.1415)*rG1112-sin(1.0*3.1415)*iG1112;
temp2 = sin(1.0*3.1415)*iG1112+cos(1.0*3.1415)*rG1112;
rG1112 = temp1;
iG1112 = temp2;

temp1 = cos(2.0*1.0*3.1415)*rG1122-sin(2.0*1.0*3.1415)*iG1122;
temp2 = sin(2.0*1.0*3.1415)*iG1122+cos(2.0*1.0*3.1415)*rG1122;
rG1122 = temp1;
iG1122 = temp2;
}




 /* THIS VERSION IS THE NEW (AND POSSIBLY CORRECT) VERSION OF THE FOCK-SPACE HAMILTONIAN */

   /* REAL PART */
   for(i=0;i<=2*np;i+=2){
     mu = (double)i/2;
     H[i][i] = (E11*N-E0+g/2.0*(G1111*(mu*mu+(N-mu)*(N-mu)-N) + G1221*(4.0*mu*(N-mu))));
     H[i+1][i+1] = H[i][i];
   }
   for(i=0;i<=2*np-2;i+=2){
     mu = (double)i/2;
     H[i][i+2] = (rE12+g*rG1112*(N-1.0))*sqrt((mu+1.0)*(N-mu));
     H[i+1][i+3] = H[i][i+2];

     H[i+2][i] = H[i][i+2];
     H[i+3][i+1] = H[i+1][i+3];
   }
   for(i=0;i<=2*np-4;i+=2){
     mu = (double)i/2;
     H[i][i+4] = g*rG1122/2.0*sqrt((N-mu-1.0)*(N-mu)*(mu+1.0)*(mu+2.0));
     H[i+1][i+5] = H[i][i+4];

     H[i+4][i] = H[i][i+4];
     H[i+5][i+1] = H[i+1][i+5];
   }   

   /* IMAGINARY PART */
   for(i=0;i<=2*np;i+=2){
     mu = (double)i/2;
     H[i+1][i] = 0.000;
     H[i][i+1] = 0.000;
   }
   for(i=0;i<=2*np-2;i+=2){
     mu = (double)i/2;
     H[i+1][i+2] = (iE12+g*iG1112*(2*mu-N+1.0))*sqrt((mu+1.0)*(N-mu));
     H[i][i+3] = -1.0*H[i+1][i+2];

     H[i+2][i+1] = 1.0*H[i+1][i+2];  
     H[i+3][i] = 1.0*H[i][i+3];
   }
   for(i=0;i<=2*np-4;i+=2){
     mu = (double)i/2;
     H[i+1][i+4] = g*iG1122/2.0*sqrt((N-mu-1.0)*(N-mu)*(mu+1.0)*(mu+2.0));
     H[i][i+5] = -1.0*H[i+1][i+4];

     H[i+4][i+1] = 1.0*H[i+1][i+4]; 
     H[i+5][i] = 1.0*H[i][i+5];
   }   

//  fprintf(stdout, "\n H[N/2][N/2] = %f \n", H[np][np]);

   /* 
   int j;
  for(i=0;i<=2*np+1;i++)
    {
    for(j=0;j<=2*np+1;j++)
      {
	fprintf(stdout, "   %e   ", H[i][j]);
      }
    fprintf(stdout, " \n ");
    }
   */

   // printf("CI basis hamiltonian filled\n");

}

void fillEngys(double engys[], double y[], double t, int nv)  // compute E12, E11 and pass them back in the engys[] vector
{

  extern double g, g12, L, N, scale_factor, *ke;
  double E11, rE12, iE12, dx, xx, mu, DFTn;
  int i;
  double *a_r, *a_i, *aTr, *aTi;
  double hbar = 63.5077;

  double extpot(double t, double xx);
  void sinft(double [], int);

  dx = L/((double)nv);
  DFTn = 2.0/((double)nv);
  a_r = vector((long)1,(long)(nv));
  a_i = vector((long)1,(long)(nv));

  // entered in x-space 

	for(i=2; i<=(2*nv); i+=2)
	{
	  a_i[i/2] = y[i];
	  a_r[i/2] = y[i-1];
	}

	E11 = 0.00;
	rE12 = 0.00;
	iE12 = 0.00;	
	xx = 0.00;
	for(i=1; i<=nv; i++)
	{
	  mu = (extpot(t, xx)/scale_factor);  
	  E11  += mu*(a_r[i]*a_r[i]+a_i[i]*a_i[i]);
	  rE12 +=  mu*(a_r[i]*a_r[nv+1-i]+a_i[i]*a_i[nv+1-i]);
	  xx += dx;
	}

	sinft(a_r, nv);
	sinft(a_i, nv);
	for(i=1;i<=nv;i++)
	{
	  a_r[i] *= ke[i]*DFTn;
	  a_i[i] *= ke[i]*DFTn;
	}	
	sinft(a_r, nv);
	sinft(a_i, nv);

	for(i=1;i<=nv;i++)
	{
	  E11 += (a_r[i]*y[2*i-1]+a_i[i]*y[2*i]);
	  rE12 += (y[2*i-1]*a_r[nv-i+1]+y[2*i]*a_i[nv+1-i]);
        }
	E11 *= dx*hbar; rE12 *= dx*hbar; iE12 *= dx*hbar;
	engys[0] = E11;
	engys[1] = rE12;
	engys[2] = iE12;

	free_vector(a_r,(long)1,(long)(nv));
	free_vector(a_i,(long)1,(long)(nv));	

}


