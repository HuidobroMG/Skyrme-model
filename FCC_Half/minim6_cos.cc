#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <time.h>

#include <cmath>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_rng.h>

#include "fields_cos.h"

#include "lineparser.h"

#define RANDOM ;

#ifdef WITHGPERFTOOLS
#include <gperftools/profiler.h>
#endif


// parameters for the qag subroutine
double epsilon = 1.e-7;
double epsabs = 1.e-3;
int key = 6;

const double bigRes = 1.e10; //1.e30 Initially


// PHYSICAL PARAMETERS
double fpi = 122; // Pion decay constant (MeV)
double e = 4.54; // Skyrme constant
double mpi = 138.0; // Pion mass (MeV)
double lambda2 = 1.0; // Sextic term coupling constant (MeV*fm**3)
double hbarc = 197.3269804064626; // hbar*c (MeV*fm)


void error(std::string pname = " " )
{
    std::cout << "Usage : " << std::endl;
    std::cout << pname << " -LLength -cCODE "
	      << std::endl;
    exit(1);
}


void ReadFile(std::string read, double * ypar, int ndim,int iBeta,int ndim2,
	      int iBeta2)
{
  
  std::cout << " Reading coefficients from file " << read << std::endl;
  std::ifstream cfile(read.c_str());
  std::vector<double> apar;
  char line[256];

  cfile.get(line,256);
  while (!cfile.eof())
    {
      int i; double a;
      cfile >> i >> a >> std::ws;
      apar.push_back(a);
    }
  cfile.close();
  if (int(apar.size()) != ndim2 )
    {
      std::cout << " Wrong coefficient file. Dimension: " << apar.size()
		<< "   expected " << ndim2 << std::endl;
      exit(1);
    }
  for (int i = 0; i< iBeta2; i++)
    {
      ypar[i] = apar[i];
    }
  for (int i = iBeta2; i< iBeta; i++)
    {
      ypar[i] = 0.0;
    }
  for (int i= iBeta2 ; i < ndim2; i++)
    {
      ypar[i-iBeta2+iBeta] = apar[i];
    }
  for (int i = ndim2-iBeta2+iBeta; i< ndim; i++)
    {
      ypar[i] = 0.0;
    }

  for (int i=0;i<ndim; i++)
    {
      std::cout << i << " " << ypar[i] << std::endl; 
    }
}

struct oparams
{
  int ndim; 
  double L;
  double a0;
  double a2;
  double a4;
  double a6;
  int c0,c2,c4,c6;
};


// Levi-Civita symbol function
double  levi3(int i, int j, int k)
{
    if((i==j)||(i==k)||(j==k))
      {
	return 0;
      }
    else
      {
	double z = static_cast<double>((j-i)*(k-i)*(k-j))/2.;
	return z;
      }
}


double  levi(int i, int j,int k, int l)
{
    if((i==j)||(i==k)||(i==l)||(j==k)||(j==l)||(k==l))
      {
	return 0;
      }
    else
      {
	double z = static_cast<double>((i-j)*(i-k)*(i-l)*(j-k)*(j-l)*(k-l))/12.;
	return z;
      }
}

// ======================================================================

// ENERGY DENSITY
double E(double *ypar,double x,double y,double z,oparams opar)
{

  std::vector<double> sig(4);
  
  double ycos[9];
  double ysin[9];
  double L = opar.L;

  ycos[0] = cos(M_PI*x/L);
  ycos[1] = cos(M_PI*y/L);
  ycos[2] = cos(M_PI*z/L);

  ycos[3] = cos(2.0*M_PI*x/L);
  ycos[4] = cos(2.0*M_PI*y/L);
  ycos[5] = cos(2.0*M_PI*z/L);
  
  ycos[6] = cos(3.0*M_PI*x/L);
  ycos[7] = cos(3.0*M_PI*y/L);
  ycos[8] = cos(3.0*M_PI*z/L);


  ysin[0] = sin(M_PI*x/L);
  ysin[1] = sin(M_PI*y/L);
  ysin[2] = sin(M_PI*z/L);

  ysin[3] = sin(2.0*M_PI*x/L);
  ysin[4] = sin(2.0*M_PI*y/L);
  ysin[5] = sin(2.0*M_PI*z/L);
  
  ysin[6] = sin(3.0*M_PI*x/L);
  ysin[7] = sin(3.0*M_PI*y/L);
  ysin[8] = sin(3.0*M_PI*z/L);
 
  
  double sigma = Sigma(ypar,ycos,ysin,L);
  double pi1   = Pi1(ypar,ycos,ysin,L);
  double pi2   = Pi2(ypar,ycos,ysin,L);
  double pi3   = Pi3(ypar,ycos,ysin,L);

  sig[0] = sigma;
  sig[1] = pi1;
  sig[2] = pi2;
  sig[3] = pi3;

  double norma2 = sigma*sigma+pi1*pi1+pi2*pi2+pi3*pi3;
  double norma4 = norma2*norma2;
  double norma = sqrt(norma2);
  
  std::vector< std::vector<double> > dsig(4);
  
  double sigma_1 = Sigma_1(ypar,ycos,ysin,L);
  double sigma_2 = Sigma_2(ypar,ycos,ysin,L);
  double sigma_3 = Sigma_3(ypar,ycos,ysin,L);
  dsig[0].push_back(sigma_1);
  dsig[0].push_back(sigma_2);
  dsig[0].push_back(sigma_3);
  
  double pi1_1 = Pi1_1(ypar,ycos,ysin,L);
  double pi1_2 = Pi1_2(ypar,ycos,ysin,L);
  double pi1_3 = Pi1_3(ypar,ycos,ysin,L);
  dsig[1].push_back(pi1_1);
  dsig[1].push_back(pi1_2);
  dsig[1].push_back(pi1_3);

  double pi2_1 = Pi2_1(ypar,ycos,ysin,L);
  double pi2_2 = Pi2_2(ypar,ycos,ysin,L);
  double pi2_3 = Pi2_3(ypar,ycos,ysin,L);
  dsig[2].push_back(pi2_1);
  dsig[2].push_back(pi2_2);
  dsig[2].push_back(pi2_3);

  double pi3_1 = Pi3_1(ypar,ycos,ysin,L);
  double pi3_2 = Pi3_2(ypar,ycos,ysin,L);
  double pi3_3 = Pi3_3(ypar,ycos,ysin,L);
  dsig[3].push_back(pi3_1);
  dsig[3].push_back(pi3_2);
  dsig[3].push_back(pi3_3);

  std::vector< std::vector<double> > dn(4);

  std::vector<double > jvec(3);
  for (int j=0;j<3;j++)
    {
      jvec[j] = 0;
      for (int i=0;i<int(sig.size());i++)
	{
	  jvec[j] = jvec[j] + sig[i]*dsig[i][j];
	}
      jvec[j] = jvec[j]/norma2;
    }

  for (int i=0;i<int(sig.size());i++)
    {
      for (int j=0;j<3; j++)
	{
	  dn[i].push_back((dsig[i][j]-sig[i]*jvec[j])/norma);
	}
    }
  
  double E2p = 0;
  double E4p = 0;
  double E6 = 0;
  for (int j = 0; j<3;j++)
    {
      for (int i=0; i<int(sig.size());i++)
	{
	  E2p = E2p + dn[i][j]*dn[i][j];
	  for (int k=0;k<3;k++)
	    {
	      for (int ip = 0;ip<int(sig.size());ip++)
		{
		  E4p = E4p +  dn[i][j]*dn[i][k]*dn[ip][j]*dn[ip][k];
		}
	    }
	}
    }

   // The sextic term
   if (opar.c6 != 0)
   {
     for (int i = 0; i<int(sig.size()); i++)
    {
      for (int j = 0; j<int(sig.size()); j++)
     {
       for (int k = 0; k<int(sig.size()); k++)
      {
	for (int g = 0; g<int(sig.size()); g++)
       {
	 //	 std::cout << "Levi " << levi(i,j,k,g) << std::endl;
	 E6 = E6 + levi(i,j,k,g)
	   *sig[i]/norma*dn[j][0]*dn[k][1]*dn[g][2];
       }
      }
     }
    }
   }
   
  double ncons = 1.0/(24.0*M_PI*M_PI);
  
  double E0 = 1.0-sigma/norma; // Pion mass term
  double E4 = E2p*E2p-E4p; // Quartic term

  double coe6 = opar.c6*2*lambda2*pow(e*e*fpi,2)/pow(hbarc,3);
  double coe0 = opar.c0*2*mpi*mpi/(e*e*fpi*fpi);
  double coe02 = 0.0*opar.c0;
  
  // Energy Density
  double Ep = ncons*(opar.c2*E2p + 2.0*opar.c4*E4 + coe6*E6*E6 + coe0*E0 + coe02*(1-sig[0]/norma)*(1-sig[0]/norma));

  return Ep;
}

//=======================================================================

// BARYON DENSITY
double B(double *ypar,double x,double y,double z,oparams opar)
{

  std::vector<double> sig(4);
  
  double ycos[9];
  double ysin[9];
  double L = opar.L;

  ycos[0] = cos(M_PI*x/L);
  ycos[1] = cos(M_PI*y/L);
  ycos[2] = cos(M_PI*z/L);

  ycos[3] = cos(2.0*M_PI*x/L);
  ycos[4] = cos(2.0*M_PI*y/L);
  ycos[5] = cos(2.0*M_PI*z/L);
  
  ycos[6] = cos(3.0*M_PI*x/L);
  ycos[7] = cos(3.0*M_PI*y/L);
  ycos[8] = cos(3.0*M_PI*z/L);


  ysin[0] = sin(M_PI*x/L);
  ysin[1] = sin(M_PI*y/L);
  ysin[2] = sin(M_PI*z/L);

  ysin[3] = sin(2.0*M_PI*x/L);
  ysin[4] = sin(2.0*M_PI*y/L);
  ysin[5] = sin(2.0*M_PI*z/L);
  
  ysin[6] = sin(3.0*M_PI*x/L);
  ysin[7] = sin(3.0*M_PI*y/L);
  ysin[8] = sin(3.0*M_PI*z/L);
 
  
  double sigma = Sigma(ypar,ycos,ysin,L);
  double pi1   = Pi1(ypar,ycos,ysin,L);
  double pi2   = Pi2(ypar,ycos,ysin,L);
  double pi3   = Pi3(ypar,ycos,ysin,L);

  sig[0] = sigma;
  sig[1] = pi1;
  sig[2] = pi2;
  sig[3] = pi3;

  double norma2 = sigma*sigma+pi1*pi1+pi2*pi2+pi3*pi3;
  double norma4 = norma2*norma2;
  double norma = sqrt(norma2);
  
  std::vector< std::vector<double> > dsig(4);
  
  double sigma_1 = Sigma_1(ypar,ycos,ysin,L);
  double sigma_2 = Sigma_2(ypar,ycos,ysin,L);
  double sigma_3 = Sigma_3(ypar,ycos,ysin,L);
  dsig[0].push_back(sigma_1);
  dsig[0].push_back(sigma_2);
  dsig[0].push_back(sigma_3);
  
  double pi1_1 = Pi1_1(ypar,ycos,ysin,L);
  double pi1_2 = Pi1_2(ypar,ycos,ysin,L);
  double pi1_3 = Pi1_3(ypar,ycos,ysin,L);
  dsig[1].push_back(pi1_1);
  dsig[1].push_back(pi1_2);
  dsig[1].push_back(pi1_3);

  double pi2_1 = Pi2_1(ypar,ycos,ysin,L);
  double pi2_2 = Pi2_2(ypar,ycos,ysin,L);
  double pi2_3 = Pi2_3(ypar,ycos,ysin,L);
  dsig[2].push_back(pi2_1);
  dsig[2].push_back(pi2_2);
  dsig[2].push_back(pi2_3);

  double pi3_1 = Pi3_1(ypar,ycos,ysin,L);
  double pi3_2 = Pi3_2(ypar,ycos,ysin,L);
  double pi3_3 = Pi3_3(ypar,ycos,ysin,L);
  dsig[3].push_back(pi3_1);
  dsig[3].push_back(pi3_2);
  dsig[3].push_back(pi3_3);

  std::vector< std::vector<double> > dn(4);

  std::vector<double > jvec(3);
  for (int j=0;j<3;j++)
    {
      jvec[j] = 0;
      for (int i=0;i<int(sig.size());i++)
	{
	  jvec[j] = jvec[j] + sig[i]*dsig[i][j];
	}
      jvec[j] = jvec[j]/norma2;
    }

  for (int i=0;i<int(sig.size());i++)
    {
      for (int j=0;j<3; j++)
	{
	  dn[i].push_back((dsig[i][j]-sig[i]*jvec[j])/norma);
	}
    }

   double B6 = 0.0;

   for (int i = 0; i<int(sig.size()); i++)
   {
    for (int j = 0; j<int(sig.size()); j++)
    {
     for (int k = 0; k<int(sig.size()); k++)
     {
      for (int g = 0; g<int(sig.size()); g++)
      {
	//	 std::cout << "Levi " << levi(i,j,k,g) << std::endl;
	B6 = B6 + levi(i,j,k,g)
	  *sig[i]/norma*dn[j][0]*dn[k][1]*dn[g][2];
      }
     }
    }
   }

   
  double ncons = 1.0/(2.0*M_PI*M_PI);

  double Bn = -ncons*B6;
  
  return Bn;
}

//=============================================================================

// ISOSPIN CONTRIBUTION
double Iso(double *ypar,double x,double y,double z,oparams opar)
{

  std::vector<double> sig(4);
  
  double ycos[9];
  double ysin[9];
  double L = opar.L;

  ycos[0] = cos(M_PI*x/L);
  ycos[1] = cos(M_PI*y/L);
  ycos[2] = cos(M_PI*z/L);

  ycos[3] = cos(2.0*M_PI*x/L);
  ycos[4] = cos(2.0*M_PI*y/L);
  ycos[5] = cos(2.0*M_PI*z/L);
  
  ycos[6] = cos(3.0*M_PI*x/L);
  ycos[7] = cos(3.0*M_PI*y/L);
  ycos[8] = cos(3.0*M_PI*z/L);


  ysin[0] = sin(M_PI*x/L);
  ysin[1] = sin(M_PI*y/L);
  ysin[2] = sin(M_PI*z/L);

  ysin[3] = sin(2.0*M_PI*x/L);
  ysin[4] = sin(2.0*M_PI*y/L);
  ysin[5] = sin(2.0*M_PI*z/L);
  
  ysin[6] = sin(3.0*M_PI*x/L);
  ysin[7] = sin(3.0*M_PI*y/L);
  ysin[8] = sin(3.0*M_PI*z/L);
  
  double sigma = Sigma(ypar,ycos,ysin,L);
  double pi1   = Pi1(ypar,ycos,ysin,L);
  double pi2   = Pi2(ypar,ycos,ysin,L);
  double pi3   = Pi3(ypar,ycos,ysin,L);

  sig[0] = sigma;
  sig[1] = pi1;
  sig[2] = pi2;
  sig[3] = pi3;

  double norma2 = sigma*sigma+pi1*pi1+pi2*pi2+pi3*pi3;
  double norma4 = norma2*norma2;
  double norma = sqrt(norma2);
  
  std::vector< std::vector<double> > dsig(4);
  
  double sigma_1 = Sigma_1(ypar,ycos,ysin,L);
  double sigma_2 = Sigma_2(ypar,ycos,ysin,L);
  double sigma_3 = Sigma_3(ypar,ycos,ysin,L);
  dsig[0].push_back(sigma_1);
  dsig[0].push_back(sigma_2);
  dsig[0].push_back(sigma_3);
  
  double pi1_1 = Pi1_1(ypar,ycos,ysin,L);
  double pi1_2 = Pi1_2(ypar,ycos,ysin,L);
  double pi1_3 = Pi1_3(ypar,ycos,ysin,L);
  dsig[1].push_back(pi1_1);
  dsig[1].push_back(pi1_2);
  dsig[1].push_back(pi1_3);

  double pi2_1 = Pi2_1(ypar,ycos,ysin,L);
  double pi2_2 = Pi2_2(ypar,ycos,ysin,L);
  double pi2_3 = Pi2_3(ypar,ycos,ysin,L);
  dsig[2].push_back(pi2_1);
  dsig[2].push_back(pi2_2);
  dsig[2].push_back(pi2_3);

  double pi3_1 = Pi3_1(ypar,ycos,ysin,L);
  double pi3_2 = Pi3_2(ypar,ycos,ysin,L);
  double pi3_3 = Pi3_3(ypar,ycos,ysin,L);
  dsig[3].push_back(pi3_1);
  dsig[3].push_back(pi3_2);
  dsig[3].push_back(pi3_3);

  std::vector< std::vector<double> > dn(4);

  std::vector<double > jvec(3);
  for (int j=0;j<3;j++)
    {
      jvec[j] = 0;
      for (int i=0;i<int(sig.size());i++)
	{
	  jvec[j] = jvec[j] + sig[i]*dsig[i][j];
	}
      jvec[j] = jvec[j]/norma2;
    }

  for (int i=0;i<int(sig.size());i++)
    {
      for (int j=0;j<3; j++)
	{
	  dn[i].push_back((dsig[i][j]-sig[i]*jvec[j])/norma);
	}
    }


   double Lambda2 = 2*(sig[1]*sig[1] + sig[2]*sig[2])/norma2;
   
   double A1 = 0;
   for (int i = 0; i < 3; i++) // Field index
   {
     A1 += dn[0][i]*dn[0][i]*(1-sig[3]*sig[3]/norma2) + dn[3][i]*dn[3][i]*(1-sig[0]*sig[0]/norma2) + 2*sig[0]*sig[3]*dn[0][i]*dn[3][i]/norma2;
   }
   double Lambda4 = 8*A1;
   double Lambda6 = 0;
   for (int i = 0; i < 3; i++) // Derivative index
  {
    for (int j = 0; j < 3; j++)
    {
      Lambda6 += (dn[0][i]*dn[3][j] - dn[0][j]*dn[3][i])*(dn[0][i]*dn[3][j] - dn[0][j]*dn[3][i]);
    }
  }
  
  double ncons = 1.0/(24.0*M_PI*M_PI);
  
  double coe6 = opar.c6*2*lambda2*fpi*fpi*pow(e,4)/pow(hbarc, 3);

  double Lambda = ncons*(Lambda2 + Lambda4 + coe6*Lambda6);

  return Lambda;
}

// J0_3 (For the charge density)
double J03(double *ypar,double x,double y,double z,oparams opar)
{

  std::vector<double> sig(4);
  
  double ycos[9];
  double ysin[9];
  double L = opar.L;

  ycos[0] = cos(M_PI*x/L);
  ycos[1] = cos(M_PI*y/L);
  ycos[2] = cos(M_PI*z/L);

  ycos[3] = cos(2.0*M_PI*x/L);
  ycos[4] = cos(2.0*M_PI*y/L);
  ycos[5] = cos(2.0*M_PI*z/L);
  
  ycos[6] = cos(3.0*M_PI*x/L);
  ycos[7] = cos(3.0*M_PI*y/L);
  ycos[8] = cos(3.0*M_PI*z/L);


  ysin[0] = sin(M_PI*x/L);
  ysin[1] = sin(M_PI*y/L);
  ysin[2] = sin(M_PI*z/L);

  ysin[3] = sin(2.0*M_PI*x/L);
  ysin[4] = sin(2.0*M_PI*y/L);
  ysin[5] = sin(2.0*M_PI*z/L);
  
  ysin[6] = sin(3.0*M_PI*x/L);
  ysin[7] = sin(3.0*M_PI*y/L);
  ysin[8] = sin(3.0*M_PI*z/L);
  
  double sigma = Sigma(ypar,ycos,ysin,L);
  double pi1   = Pi1(ypar,ycos,ysin,L);
  double pi2   = Pi2(ypar,ycos,ysin,L);
  double pi3   = Pi3(ypar,ycos,ysin,L);

  sig[0] = sigma;
  sig[1] = pi1;
  sig[2] = pi2;
  sig[3] = pi3;

  double norma2 = sigma*sigma+pi1*pi1+pi2*pi2+pi3*pi3;
  double norma4 = norma2*norma2;
  double norma = sqrt(norma2);
  
  std::vector< std::vector<double> > dsig(4);
  
  double sigma_1 = Sigma_1(ypar,ycos,ysin,L);
  double sigma_2 = Sigma_2(ypar,ycos,ysin,L);
  double sigma_3 = Sigma_3(ypar,ycos,ysin,L);
  dsig[0].push_back(sigma_1);
  dsig[0].push_back(sigma_2);
  dsig[0].push_back(sigma_3);
  
  double pi1_1 = Pi1_1(ypar,ycos,ysin,L);
  double pi1_2 = Pi1_2(ypar,ycos,ysin,L);
  double pi1_3 = Pi1_3(ypar,ycos,ysin,L);
  dsig[1].push_back(pi1_1);
  dsig[1].push_back(pi1_2);
  dsig[1].push_back(pi1_3);

  double pi2_1 = Pi2_1(ypar,ycos,ysin,L);
  double pi2_2 = Pi2_2(ypar,ycos,ysin,L);
  double pi2_3 = Pi2_3(ypar,ycos,ysin,L);
  dsig[2].push_back(pi2_1);
  dsig[2].push_back(pi2_2);
  dsig[2].push_back(pi2_3);

  double pi3_1 = Pi3_1(ypar,ycos,ysin,L);
  double pi3_2 = Pi3_2(ypar,ycos,ysin,L);
  double pi3_3 = Pi3_3(ypar,ycos,ysin,L);
  dsig[3].push_back(pi3_1);
  dsig[3].push_back(pi3_2);
  dsig[3].push_back(pi3_3);

  std::vector< std::vector<double> > dn(4);

  std::vector<double > jvec(3);
  for (int j=0;j<3;j++)
    {
      jvec[j] = 0;
      for (int i=0;i<int(sig.size());i++)
	{
	  jvec[j] = jvec[j] + sig[i]*dsig[i][j];
	}
      jvec[j] = jvec[j]/norma2;
    }

  for (int i=0;i<int(sig.size());i++)
    {
      for (int j=0;j<3; j++)
	{
	  dn[i].push_back((dsig[i][j]-sig[i]*jvec[j])/norma);
	}
    }
  
  // Expected value matrix of < Ri Kj > in the ground state
  double RiKj[4][4] = {{0, 0, 0, 0}, {0, 2.0/3.0, 1.0/3.0, 1.0/3.0}, {0, 1.0/3.0, 2.0/3.0, -1.0/3.0}, {0, 1.0/3.0, -1.0/3.0, 2.0/3.0}};
  
  double delta[4][4] = {{0, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
  
  double ncons = 1/(24.0*M_PI*M_PI);
  
  // Quartic term
  double l4RK_1 = 0;
  double l4RK_2 = 0;
  double l4RK_3 = 0;
  double l4RK_4 = 0;
  double l4RK_5 = 0;
  double l4RK_6 = 0;
  double l4RK_7 = 0;
  double l4RK_8 = 0;
  double l4RK_9 = 0;
  for (int i = 1; i < int(sig.size()); i++)
  {
    for (int j = 1; j < int(sig.size()); j++)
    {
      l4RK_1 += (sig[i]*sig[j] - (sig[1]*sig[1] + sig[2]*sig[2] + sig[3]*sig[3])*delta[i][j])/norma2*RiKj[i][j];
      l4RK_7 += sig[i]*sig[j]/norma2*RiKj[i][j];
      
      for (int k = 0; k < 3; k++)
      {
        l4RK_8 += dn[i][k]*dn[j][k]*RiKj[i][j];

        for (int a = 1; a < int(sig.size()); a++)
        {
          l4RK_9 += (sig[a]*dn[a][k])*(sig[i]*dn[j][k] + sig[j]*dn[i][k])/norma2*RiKj[i][j];
          
          for (int b = 1; b < int(sig.size()); b++)
          {
            l4RK_4 += levi3(i,a,b)*sig[a]*dn[b][k]/norma*(sig[1]*sig[1] + sig[2]*sig[2] + sig[3]*sig[3])*dn[j][k]/norma2*RiKj[i][j];
            l4RK_6 += levi3(j,a,b)*sig[a]*dn[b][k]/norma*(sig[1]*sig[1] + sig[2]*sig[2] + sig[3]*sig[3])*dn[i][k]/norma2*RiKj[i][j];
            
            for (int c = 1; c < int(sig.size()); c++)
            {
              l4RK_3 += levi3(i,a,b)*sig[a]*dn[b][k]*sig[j]/norma2*sig[c]*dn[c][k]/norma*RiKj[i][j];
              l4RK_5 += levi3(j,a,b)*sig[a]*dn[b][k]*sig[i]/norma2*sig[c]*dn[c][k]/norma*RiKj[i][j];
              
              for (int d = 1; d < int(sig.size()); d++)
              {
                l4RK_2 += levi3(i,a,b)*sig[a]*dn[b][k]*levi3(j,c,d)*sig[c]*dn[d][k]/norma2*RiKj[i][j];
              }
            }
          }
        }
      }
    }
  }
    
  double A1 = 0;
  double A2 = 0;
  double A3 = 0;
  double A4 = 0;
  double A5 = 0;
  for (int i = 1; i<int(sig.size()); i++) // Field index
  {
    A1 += sig[i]*sig[i]/norma2;
    for (int k = 0; k<3; k++)  // Derivative index
      {
        A2 += dn[i][k]*dn[i][k];
        A3 += sig[i]*sig[i]*dn[0][k]*dn[0][k]/norma2;
        A4 += sig[0]*dn[0][k]*sig[i]*dn[i][k]/norma2;
        for (int j = 1; j<int(sig.size()); j++)
        {
          A5 += sig[i]*dn[i][k]*sig[j]*dn[j][k]/norma2;
        }
      }
  }
  
  double l4RK = 8*l4RK_1*(sig[0]*sig[0]/norma2*A2 + A3 - 2*A4 + A1*A2 - A5) + 8*l4RK_2 + 16*sig[0]/norma*(l4RK_3 - l4RK_4 + l4RK_5 - l4RK_6) + 32.0*sig[0]*sig[0]/norma2*(A5*l4RK_7 + A1*A1*l4RK_8 - A1*l4RK_9);
  
  // Quadratic term
  double l2RK = 2*l4RK_1;
  
  return ncons*(-l2RK - l4RK);
  
}

//===============================================================================

// Kaon condensate Potential
double Vk_Cond(double *ypar,double x,double y,double z,oparams opar)
{

  std::vector<double> sig(4);
  
  double ycos[9];
  double ysin[9];
  double L = opar.L;

  ycos[0] = cos(M_PI*x/L);
  ycos[1] = cos(M_PI*y/L);
  ycos[2] = cos(M_PI*z/L);

  ycos[3] = cos(2.0*M_PI*x/L);
  ycos[4] = cos(2.0*M_PI*y/L);
  ycos[5] = cos(2.0*M_PI*z/L);
  
  ycos[6] = cos(3.0*M_PI*x/L);
  ycos[7] = cos(3.0*M_PI*y/L);
  ycos[8] = cos(3.0*M_PI*z/L);


  ysin[0] = sin(M_PI*x/L);
  ysin[1] = sin(M_PI*y/L);
  ysin[2] = sin(M_PI*z/L);

  ysin[3] = sin(2.0*M_PI*x/L);
  ysin[4] = sin(2.0*M_PI*y/L);
  ysin[5] = sin(2.0*M_PI*z/L);
  
  ysin[6] = sin(3.0*M_PI*x/L);
  ysin[7] = sin(3.0*M_PI*y/L);
  ysin[8] = sin(3.0*M_PI*z/L);
 
  
  double sigma = Sigma(ypar,ycos,ysin,L);
  double pi1   = Pi1(ypar,ycos,ysin,L);
  double pi2   = Pi2(ypar,ycos,ysin,L);
  double pi3   = Pi3(ypar,ycos,ysin,L);

  sig[0] = sigma;
  sig[1] = pi1;
  sig[2] = pi2;
  sig[3] = pi3;

  double norma2 = sigma*sigma+pi1*pi1+pi2*pi2+pi3*pi3;
  double norma4 = norma2*norma2;
  double norma = sqrt(norma2);
  
  std::vector< std::vector<double> > dsig(4);
  
  double sigma_1 = Sigma_1(ypar,ycos,ysin,L);
  double sigma_2 = Sigma_2(ypar,ycos,ysin,L);
  double sigma_3 = Sigma_3(ypar,ycos,ysin,L);
  dsig[0].push_back(sigma_1);
  dsig[0].push_back(sigma_2);
  dsig[0].push_back(sigma_3);
  
  double pi1_1 = Pi1_1(ypar,ycos,ysin,L);
  double pi1_2 = Pi1_2(ypar,ycos,ysin,L);
  double pi1_3 = Pi1_3(ypar,ycos,ysin,L);
  dsig[1].push_back(pi1_1);
  dsig[1].push_back(pi1_2);
  dsig[1].push_back(pi1_3);

  double pi2_1 = Pi2_1(ypar,ycos,ysin,L);
  double pi2_2 = Pi2_2(ypar,ycos,ysin,L);
  double pi2_3 = Pi2_3(ypar,ycos,ysin,L);
  dsig[2].push_back(pi2_1);
  dsig[2].push_back(pi2_2);
  dsig[2].push_back(pi2_3);

  double pi3_1 = Pi3_1(ypar,ycos,ysin,L);
  double pi3_2 = Pi3_2(ypar,ycos,ysin,L);
  double pi3_3 = Pi3_3(ypar,ycos,ysin,L);
  dsig[3].push_back(pi3_1);
  dsig[3].push_back(pi3_2);
  dsig[3].push_back(pi3_3);

  std::vector< std::vector<double> > dn(4);

  std::vector<double > jvec(3);
  for (int j=0;j<3;j++)
    {
      jvec[j] = 0;
      for (int i=0;i<int(sig.size());i++)
	{
	  jvec[j] = jvec[j] + sig[i]*dsig[i][j];
	}
      jvec[j] = jvec[j]/norma2;
    }

  for (int i=0;i<int(sig.size());i++)
    {
      for (int j=0;j<3; j++)
	{
	  dn[i].push_back((dsig[i][j]-sig[i]*jvec[j])/norma);
	}
    }

  double A1 = 1.0 + (sig[0]*sig[0] + sig[3]*sig[3])/norma2; // 1 + sigma**2 + pi3**2
  double A2 = 0; // (1 + sigma) * d_i n**2
  double A3 = 0; // d_i sigma**2 * (1 - pi3**2) + d_i pi3**2 * (1 - sigma**2) + 2 * sigma * d_i sigma * pi3 * di_ pi3 
  double A4 = 0; // Sextic term part. A_ij*A_ij
  
  for (int k = 0; k<3; k++){
       A2 += (1.0 + sig[0]/norma)*(dn[0][k]*dn[0][k] + dn[1][k]*dn[1][k] + dn[2][k]*dn[2][k] + dn[3][k]*dn[3][k]);
       A3 += dn[0][k]*dn[0][k]*(1-sig[3]*sig[3]/norma2) + dn[3][k]*dn[3][k]*(1-sig[0]*sig[0]/norma2) + 2*dn[0][k]*dn[3][k]*sig[0]*sig[3]/norma2;
  }

   for (int i = 0; i < 3; i++) // Derivative index
  {
    for (int j = 0; j < 3; j++)
    {
      A4 += (dn[0][i]*dn[3][j] - dn[0][j]*dn[3][i])*(dn[0][i]*dn[3][j] - dn[0][j]*dn[3][i]);
    }
  }
  
  
  return opar.c0*A1 + opar.c2*A2 + opar.c4*A3 + opar.c6*A4;
}


// Isospin inertia tensor of the crystal with kaons
double Lambda_K(double *ypar,double x,double y,double z,oparams opar)
{

  std::vector<double> sig(4);
  
  double ycos[9];
  double ysin[9];
  double L = opar.L;

  ycos[0] = cos(M_PI*x/L);
  ycos[1] = cos(M_PI*y/L);
  ycos[2] = cos(M_PI*z/L);

  ycos[3] = cos(2.0*M_PI*x/L);
  ycos[4] = cos(2.0*M_PI*y/L);
  ycos[5] = cos(2.0*M_PI*z/L);
  
  ycos[6] = cos(3.0*M_PI*x/L);
  ycos[7] = cos(3.0*M_PI*y/L);
  ycos[8] = cos(3.0*M_PI*z/L);


  ysin[0] = sin(M_PI*x/L);
  ysin[1] = sin(M_PI*y/L);
  ysin[2] = sin(M_PI*z/L);

  ysin[3] = sin(2.0*M_PI*x/L);
  ysin[4] = sin(2.0*M_PI*y/L);
  ysin[5] = sin(2.0*M_PI*z/L);
  
  ysin[6] = sin(3.0*M_PI*x/L);
  ysin[7] = sin(3.0*M_PI*y/L);
  ysin[8] = sin(3.0*M_PI*z/L);
 
  
  double sigma = Sigma(ypar,ycos,ysin,L);
  double pi1   = Pi1(ypar,ycos,ysin,L);
  double pi2   = Pi2(ypar,ycos,ysin,L);
  double pi3   = Pi3(ypar,ycos,ysin,L);

  sig[0] = sigma;
  sig[1] = pi1;
  sig[2] = pi2;
  sig[3] = pi3;

  double norma2 = sigma*sigma+pi1*pi1+pi2*pi2+pi3*pi3;
  double norma4 = norma2*norma2;
  double norma = sqrt(norma2);
  
  std::vector< std::vector<double> > dsig(4);
  
  double sigma_1 = Sigma_1(ypar,ycos,ysin,L);
  double sigma_2 = Sigma_2(ypar,ycos,ysin,L);
  double sigma_3 = Sigma_3(ypar,ycos,ysin,L);
  dsig[0].push_back(sigma_1);
  dsig[0].push_back(sigma_2);
  dsig[0].push_back(sigma_3);
  
  double pi1_1 = Pi1_1(ypar,ycos,ysin,L);
  double pi1_2 = Pi1_2(ypar,ycos,ysin,L);
  double pi1_3 = Pi1_3(ypar,ycos,ysin,L);
  dsig[1].push_back(pi1_1);
  dsig[1].push_back(pi1_2);
  dsig[1].push_back(pi1_3);

  double pi2_1 = Pi2_1(ypar,ycos,ysin,L);
  double pi2_2 = Pi2_2(ypar,ycos,ysin,L);
  double pi2_3 = Pi2_3(ypar,ycos,ysin,L);
  dsig[2].push_back(pi2_1);
  dsig[2].push_back(pi2_2);
  dsig[2].push_back(pi2_3);

  double pi3_1 = Pi3_1(ypar,ycos,ysin,L);
  double pi3_2 = Pi3_2(ypar,ycos,ysin,L);
  double pi3_3 = Pi3_3(ypar,ycos,ysin,L);
  dsig[3].push_back(pi3_1);
  dsig[3].push_back(pi3_2);
  dsig[3].push_back(pi3_3);

  std::vector< std::vector<double> > dn(4);

  std::vector<double > jvec(3);
  for (int j=0;j<3;j++)
    {
      jvec[j] = 0;
      for (int i=0;i<int(sig.size());i++)
	{
	  jvec[j] = jvec[j] + sig[i]*dsig[i][j];
	}
      jvec[j] = jvec[j]/norma2;
    }

  for (int i=0;i<int(sig.size());i++)
    {
      for (int j=0;j<3; j++)
	{
	  dn[i].push_back((dsig[i][j]-sig[i]*jvec[j])/norma);
	}
    }

  double L1 = (sig[1]*sig[1] + sig[2]*sig[2])/norma2; // pi_1 **2 + pi_2 **2
  double L2 = 0; // pi_1 **2 * d_i pi_2 **2 + pi_2 **2 * d_i pi_1 **2 - 2 * pi_1 * pi_2 * d_i pi_1 * d_i pi_1
  double L3 = 0; // d_i n**2 * (pi_1 **2 + pi_2 **2)
  
  for (int i = 0; i < 3; i++){
  	L2 += sig[1]*sig[1]*dn[2][i]*dn[2][i]/norma2 + sig[2]*sig[2]*dn[1][i]*dn[1][i]/norma2 - 2*sig[1]*sig[2]*dn[1][i]*dn[2][i]/norma2;
  	L3 += L1*(dn[0][i]*dn[0][i] + dn[1][i]*dn[1][i] + dn[2][i]*dn[2][i] + dn[3][i]*dn[3][i]);
  }
  
  return opar.c0*L1 + opar.c2*L2 + opar.c4*L3;
}

//===============================================================================

// Pion fluctuations

// Gamma(L) and H(L)
double GH(double *ypar,double x,double y,double z,oparams opar)
{

  std::vector<double> sig(4);
  
  double ycos[9];
  double ysin[9];
  double L = opar.L;

  ycos[0] = cos(M_PI*x/L);
  ycos[1] = cos(M_PI*y/L);
  ycos[2] = cos(M_PI*z/L);

  ycos[3] = cos(2.0*M_PI*x/L);
  ycos[4] = cos(2.0*M_PI*y/L);
  ycos[5] = cos(2.0*M_PI*z/L);
  
  ycos[6] = cos(3.0*M_PI*x/L);
  ycos[7] = cos(3.0*M_PI*y/L);
  ycos[8] = cos(3.0*M_PI*z/L);


  ysin[0] = sin(M_PI*x/L);
  ysin[1] = sin(M_PI*y/L);
  ysin[2] = sin(M_PI*z/L);

  ysin[3] = sin(2.0*M_PI*x/L);
  ysin[4] = sin(2.0*M_PI*y/L);
  ysin[5] = sin(2.0*M_PI*z/L);
  
  ysin[6] = sin(3.0*M_PI*x/L);
  ysin[7] = sin(3.0*M_PI*y/L);
  ysin[8] = sin(3.0*M_PI*z/L);
 
  
  double sigma = Sigma(ypar,ycos,ysin,L);
  double pi1   = Pi1(ypar,ycos,ysin,L);
  double pi2   = Pi2(ypar,ycos,ysin,L);
  double pi3   = Pi3(ypar,ycos,ysin,L);

  sig[0] = sigma;
  sig[1] = pi1;
  sig[2] = pi2;
  sig[3] = pi3;

  double norma2 = sigma*sigma+pi1*pi1+pi2*pi2+pi3*pi3;
  double norma4 = norma2*norma2;
  double norma = sqrt(norma2);
  
  std::vector< std::vector<double> > dsig(4);
  
  double sigma_1 = Sigma_1(ypar,ycos,ysin,L);
  double sigma_2 = Sigma_2(ypar,ycos,ysin,L);
  double sigma_3 = Sigma_3(ypar,ycos,ysin,L);
  dsig[0].push_back(sigma_1);
  dsig[0].push_back(sigma_2);
  dsig[0].push_back(sigma_3);
  
  double pi1_1 = Pi1_1(ypar,ycos,ysin,L);
  double pi1_2 = Pi1_2(ypar,ycos,ysin,L);
  double pi1_3 = Pi1_3(ypar,ycos,ysin,L);
  dsig[1].push_back(pi1_1);
  dsig[1].push_back(pi1_2);
  dsig[1].push_back(pi1_3);

  double pi2_1 = Pi2_1(ypar,ycos,ysin,L);
  double pi2_2 = Pi2_2(ypar,ycos,ysin,L);
  double pi2_3 = Pi2_3(ypar,ycos,ysin,L);
  dsig[2].push_back(pi2_1);
  dsig[2].push_back(pi2_2);
  dsig[2].push_back(pi2_3);

  double pi3_1 = Pi3_1(ypar,ycos,ysin,L);
  double pi3_2 = Pi3_2(ypar,ycos,ysin,L);
  double pi3_3 = Pi3_3(ypar,ycos,ysin,L);
  dsig[3].push_back(pi3_1);
  dsig[3].push_back(pi3_2);
  dsig[3].push_back(pi3_3);

  std::vector< std::vector<double> > dn(4);

  std::vector<double > jvec(3);
  for (int j=0;j<3;j++)
    {
      jvec[j] = 0;
      for (int i=0;i<int(sig.size());i++)
	{
	  jvec[j] = jvec[j] + sig[i]*dsig[i][j];
	}
      jvec[j] = jvec[j]/norma2;
    }

  for (int i=0;i<int(sig.size());i++)
    {
      for (int j=0;j<3; j++)
	{
	  dn[i].push_back((dsig[i][j]-sig[i]*jvec[j])/norma);
	}
    }
  
  double A1 = 0.0;
  double A2 = 0.0;
  double A3 = 0.0;
  double A4 = dn[0][0]*dn[0][0]+dn[0][1]*dn[0][1]+dn[0][2]*dn[0][2];
  for (int i=1;i<int(sig.size());i++)
   {
     A1 += sig[i]*sig[i]/norma2;
     for (int k=0;k<3;k++)
     {
       A2 += dn[i][k]*dn[i][k];
       A3 += dn[0][k]*sig[i]/norma*dn[i][k];
     }
   }
   
  double Gamma_L = -4.0/3.0*A1 + 8.0*((A4+A2)*(sig[0]*sig[0]/norma2+A1/3.0) - 1.0/3.0*(A4*A1-2.0*sig[0]/norma*A3 + sig[0]*sig[0]/norma2*A2));
  
  double H_L = -4.0/3.0*A1 + 16.0/3.0*((A4+A2)*(sig[0]*sig[0]/norma2+A1/3.0) - 1.0/3.0*(A4*A1-2.0*sig[0]/norma*A3 + sig[0]*sig[0]/norma2*A2));

  return opar.c0*Gamma_L + opar.c2*H_L;
}

//======================================================================================

struct xyparams
{
  double x;
  double y;
  double * app;
  oparams opar;
};

//==========================================================================

// E
double xyzIntegrand(double z, void * p)
{
  struct xyparams * apar = static_cast<struct xyparams *> (p);
  double x = apar->x;
  double y = apar->y;
  double * ypar = apar->app;
  oparams mopar = apar->opar;
  
  double Es = E(ypar,x,y,z,mopar);
  return Es;
}

// B
double xyzIntegrandB(double z, void * p)
{
  struct xyparams * apar = static_cast<struct xyparams *> (p);
  double x = apar->x;
  double y = apar->y;
  double * ypar = apar->app;
  oparams mopar = apar->opar;
  
  double Bs = B(ypar,x,y,z,mopar);
  return Bs;
}

// Isospin
double xyzIntegrandIso(double z, void * p)
{
  struct xyparams * apar = static_cast<struct xyparams *> (p);
  double x = apar->x;
  double y = apar->y;
  double * ypar = apar->app;
  oparams mopar = apar->opar;
  
  double Isos = Iso(ypar,x,y,z,mopar);
  return Isos;
}

// Lambdaij
double xyzIntegrandLambda(double z, void * p)
{
  struct xyparams * apar = static_cast<struct xyparams *> (p);
  double x = apar->x;
  double y = apar->y;
  double * ypar = apar->app;
  oparams mopar = apar->opar;
  
  double Lb = J03(ypar,x,y,z,mopar);
  return Lb;
}

// Kaons Condensate
double xyzIntegrandVk(double z, void * p)
{
  struct xyparams * apar = static_cast<struct xyparams *> (p);
  double x = apar->x;
  double y = apar->y;
  double * ypar = apar->app;
  oparams mopar = apar->opar;
  
  double Vk = Vk_Cond(ypar,x,y,z,mopar);
  return Vk;
}

// Isospin Inertia tensor with Kaons
double xyzIntegrandLambda_K(double z, void * p)
{
  struct xyparams * apar = static_cast<struct xyparams *> (p);
  double x = apar->x;
  double y = apar->y;
  double * ypar = apar->app;
  oparams mopar = apar->opar;
  
  double Lambdak = Lambda_K(ypar,x,y,z,mopar);
  return Lambdak;
}

// Gamma and H
double xyzIntegrandGH(double z, void * p)
{
  struct xyparams * apar = static_cast<struct xyparams *> (p);
  double x = apar->x;
  double y = apar->y;
  double * ypar = apar->app;
  oparams mopar = apar->opar;
  
  double Gammas = GH(ypar,x,y,z,mopar);
  return Gammas;
}

//==========================================================================

struct xparams
{
  double x;
  double * app;
  oparams opar;
};

//==========================================================================

// E
double xyIntegrand(double y,void* p)
{

  struct xparams * apar = static_cast<struct xparams *> (p);
  double x = apar->x;
  double *app = apar->app;
  oparams mopar = apar->opar;
  double Llim = mopar.L;
  
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  gsl_function F;
  struct xyparams  params = { x, y, app, mopar};

  F.function = &xyzIntegrand;
  F.params = &params;

  double result, error;
  gsl_set_error_handler_off();
  int ere = gsl_integration_qag (&F, -Llim, Llim, epsabs, epsilon, 1000, key, w,
				 &result,
				 &error);

  gsl_integration_workspace_free (w);
  if (ere !=0 ) result = bigRes;
  return result;
}

// B
double xyIntegrandB(double y,void* p)
{

  struct xparams * apar = static_cast<struct xparams *> (p);
  double x = apar->x;
  double *app = apar->app;
  oparams mopar = apar->opar;
  double Llim = mopar.L;
  
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  gsl_function F;
  struct xyparams  params = { x, y, app, mopar};

  F.function = &xyzIntegrandB;
  F.params = &params;

  double result, error;
  gsl_set_error_handler_off();
  int ere = gsl_integration_qag (&F, -Llim, Llim, epsabs, epsilon, 1000, key, w,
				 &result,
				 &error);

  gsl_integration_workspace_free (w);
  if (ere !=0 ) result = bigRes;
  return result;
}

// Isospin
double xyIntegrandIso(double y,void* p)
{

  struct xparams * apar = static_cast<struct xparams *> (p);
  double x = apar->x;
  double *app = apar->app;
  oparams mopar = apar->opar;
  double Llim = mopar.L;
  
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  gsl_function F;
  struct xyparams  params = { x, y, app, mopar};

  F.function = &xyzIntegrandIso;
  F.params = &params;

  double result, error;
  gsl_set_error_handler_off();
  int ere = gsl_integration_qag (&F, -Llim, Llim, epsabs, epsilon, 1000, key, w,
				 &result,
				 &error);

  gsl_integration_workspace_free (w);
  if (ere !=0 ) result = bigRes;
  return result;
}

// Lambdaij
double xyIntegrandLambda(double y,void* p)
{

  struct xparams * apar = static_cast<struct xparams *> (p);
  double x = apar->x;
  double *app = apar->app;
  oparams mopar = apar->opar;
  double Llim = mopar.L;
  
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  gsl_function F;
  struct xyparams  params = { x, y, app, mopar};

  F.function = &xyzIntegrandLambda;
  F.params = &params;

  double result, error;
  gsl_set_error_handler_off();
  int ere = gsl_integration_qag (&F, -Llim, Llim, epsabs, epsilon, 1000, key, w,
				 &result,
				 &error);

  gsl_integration_workspace_free (w);
  if (ere !=0 ) result = bigRes;
  return result;
}

// Kaons Condensate
double xyIntegrandVk(double y,void* p)
{

  struct xparams * apar = static_cast<struct xparams *> (p);
  double x = apar->x;
  double *app = apar->app;
  oparams mopar = apar->opar;
  double Llim = mopar.L;
  
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  gsl_function F;
  struct xyparams  params = { x, y, app, mopar};

  F.function = &xyzIntegrandVk;
  F.params = &params;

  double result, error;
  gsl_set_error_handler_off();
  int ere = gsl_integration_qag (&F, -Llim, Llim, epsabs, epsilon, 1000, key, w,
				 &result,
				 &error);

  gsl_integration_workspace_free (w);
  if (ere !=0 ) result = bigRes;
  return result;
}

double xyIntegrandLambda_K(double y,void* p)
{

  struct xparams * apar = static_cast<struct xparams *> (p);
  double x = apar->x;
  double *app = apar->app;
  oparams mopar = apar->opar;
  double Llim = mopar.L;
  
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  gsl_function F;
  struct xyparams  params = { x, y, app, mopar};

  F.function = &xyzIntegrandLambda_K;
  F.params = &params;

  double result, error;
  gsl_set_error_handler_off();
  int ere = gsl_integration_qag (&F, -Llim, Llim, epsabs, epsilon, 1000, key, w,
				 &result,
				 &error);

  gsl_integration_workspace_free (w);
  if (ere !=0 ) result = bigRes;
  return result;
}

// Gamma, H
double xyIntegrandGH(double y,void* p)
{

  struct xparams * apar = static_cast<struct xparams *> (p);
  double x = apar->x;
  double *app = apar->app;
  oparams mopar = apar->opar;
  double Llim = mopar.L;
  
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  gsl_function F;
  struct xyparams  params = { x, y, app, mopar};

  F.function = &xyzIntegrandGH;
  F.params = &params;

  double result, error;
  gsl_set_error_handler_off();
  int ere = gsl_integration_qag (&F, -Llim, Llim, epsabs, epsilon, 1000, key, w,
				 &result,
				 &error);

  gsl_integration_workspace_free (w);
  if (ere !=0 ) result = bigRes;
  return result;
}

//==========================================================================

struct Aparams
{
  double * app;
  oparams opar;
};

//==========================================================================

// E
double xIntegrand(double x,void* p)
{

  struct Aparams * apar = static_cast<struct Aparams *> (p);
  double *app = apar->app;
  oparams mopar = apar->opar;
  double Llim = mopar.L;
  
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  gsl_function F;
  struct xparams  params = { x, app, mopar };

  F.function = &xyIntegrand;
  F.params = &params;

  double result, error;
  gsl_set_error_handler_off();
  int ere = gsl_integration_qag (&F, -Llim, Llim, epsabs, epsilon, 1000, key, w,
				 &result,
				 &error);
  gsl_integration_workspace_free (w);

  if (ere !=0 ) result = bigRes;

  return result;
}

// B
double xIntegrandB(double x,void* p)
{

  struct Aparams * apar = static_cast<struct Aparams *> (p);
  double *app = apar->app;
  oparams mopar = apar->opar;
  double Llim = mopar.L;
  
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  gsl_function F;
  struct xparams  params = { x, app, mopar };

  F.function = &xyIntegrandB;
  F.params = &params;

  double result, error;
  gsl_set_error_handler_off();
  int ere = gsl_integration_qag (&F, -Llim, Llim, epsabs, epsilon, 1000, key, w,
				 &result,
				 &error);
  gsl_integration_workspace_free (w);

  if (ere !=0 ) result = bigRes;

  return result;
}

// Isospin
double xIntegrandIso(double x,void* p)
{

  struct Aparams * apar = static_cast<struct Aparams *> (p);
  double *app = apar->app;
  oparams mopar = apar->opar;
  double Llim = mopar.L;
  
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  gsl_function F;
  struct xparams  params = { x, app, mopar };

  F.function = &xyIntegrandIso;
  F.params = &params;

  double result, error;
  gsl_set_error_handler_off();
  int ere = gsl_integration_qag (&F, -Llim, Llim, epsabs, epsilon, 1000, key, w,
				 &result,
				 &error);
  gsl_integration_workspace_free (w);

  if (ere !=0 ) result = bigRes;

  return result;
}

// Lambdaij
double xIntegrandLambda(double x,void* p)
{

  struct Aparams * apar = static_cast<struct Aparams *> (p);
  double *app = apar->app;
  oparams mopar = apar->opar;
  double Llim = mopar.L;
  
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  gsl_function F;
  struct xparams  params = { x, app, mopar };

  F.function = &xyIntegrandLambda;
  F.params = &params;

  double result, error;
  gsl_set_error_handler_off();
  int ere = gsl_integration_qag (&F, -Llim, Llim, epsabs, epsilon, 1000, key, w,
				 &result,
				 &error);
  gsl_integration_workspace_free (w);

  if (ere !=0 ) result = bigRes;

  return result;
}

// Kaons Condensate
double xIntegrandVk(double x,void* p)
{

  struct Aparams * apar = static_cast<struct Aparams *> (p);
  double *app = apar->app;
  oparams mopar = apar->opar;
  double Llim = mopar.L;
  
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  gsl_function F;
  struct xparams  params = { x, app, mopar };

  F.function = &xyIntegrandVk;
  F.params = &params;

  double result, error;
  gsl_set_error_handler_off();
  int ere = gsl_integration_qag (&F, -Llim, Llim, epsabs, epsilon, 1000, key, w,
				 &result,
				 &error);
  gsl_integration_workspace_free (w);

  if (ere !=0 ) result = bigRes;

  return result;
}

double xIntegrandLambda_K(double x,void* p)
{

  struct Aparams * apar = static_cast<struct Aparams *> (p);
  double *app = apar->app;
  oparams mopar = apar->opar;
  double Llim = mopar.L;
  
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  gsl_function F;
  struct xparams  params = { x, app, mopar };

  F.function = &xyIntegrandLambda_K;
  F.params = &params;

  double result, error;
  gsl_set_error_handler_off();
  int ere = gsl_integration_qag (&F, -Llim, Llim, epsabs, epsilon, 1000, key, w,
				 &result,
				 &error);
  gsl_integration_workspace_free (w);

  if (ere !=0 ) result = bigRes;

  return result;
}

// Gamma
double xIntegrandGH(double x,void* p)
{

  struct Aparams * apar = static_cast<struct Aparams *> (p);
  double *app = apar->app;
  oparams mopar = apar->opar;
  double Llim = mopar.L;
  
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  gsl_function F;
  struct xparams  params = { x, app, mopar };

  F.function = &xyIntegrandGH;
  F.params = &params;

  double result, error;
  gsl_set_error_handler_off();
  int ere = gsl_integration_qag (&F, -Llim, Llim, epsabs, epsilon, 1000, key, w,
				 &result,
				 &error);
  gsl_integration_workspace_free (w);

  if (ere !=0 ) result = bigRes;

  return result;
}

//==========================================================================

// E
double Integrand(Aparams * app)
{
  
  oparams mopar = app->opar;
  double Llim = mopar.L;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  gsl_function F;

  F.function = &xIntegrand;
  F.params = app;

  double result, error;
  gsl_set_error_handler_off();
  int ere = gsl_integration_qag (&F, -Llim, Llim, epsabs, epsilon, 1000, key, w,
				 &result,
				 &error);
  gsl_integration_workspace_free (w);

  if (ere !=0 ) result = bigRes;

  return result;
}

// B
double IntegrandB(Aparams * app)
{
  
  oparams mopar = app->opar;
  double Llim = mopar.L;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  gsl_function F;

  F.function = &xIntegrandB;
  F.params = app;

  double result, error;
  gsl_set_error_handler_off();
  int ere = gsl_integration_qag (&F, -Llim, Llim, epsabs, epsilon, 1000, key, w,
				 &result,
				 &error);
  gsl_integration_workspace_free (w);

  if (ere !=0 ) result = bigRes;

  return result;
}

// Isospin
double IntegrandIso(Aparams * app)
{
  
  oparams mopar = app->opar;
  double Llim = mopar.L;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  gsl_function F;

  F.function = &xIntegrandIso;
  F.params = app;

  double result, error;
  gsl_set_error_handler_off();
  int ere = gsl_integration_qag (&F, -Llim, Llim, epsabs, epsilon, 1000, key, w,
				 &result,
				 &error);
  gsl_integration_workspace_free (w);

  if (ere !=0 ) result = bigRes;

  return result;
}

// Lambdaij
double IntegrandLambda(Aparams * app)
{
  
  oparams mopar = app->opar;
  double Llim = mopar.L;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  gsl_function F;

  F.function = &xIntegrandLambda;
  F.params = app;

  double result, error;
  gsl_set_error_handler_off();
  int ere = gsl_integration_qag (&F, -Llim, Llim, epsabs, epsilon, 1000, key, w,
				 &result,
				 &error);
  gsl_integration_workspace_free (w);

  if (ere !=0 ) result = bigRes;

  return result;
}

// Kaons Condensate
double IntegrandVk(Aparams * app)
{
  
  oparams mopar = app->opar;
  double Llim = mopar.L;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  gsl_function F;

  F.function = &xIntegrandVk;
  F.params = app;

  double result, error;
  gsl_set_error_handler_off();
  int ere = gsl_integration_qag (&F, -Llim, Llim, epsabs, epsilon, 1000, key, w,
				 &result,
				 &error);
  gsl_integration_workspace_free (w);

  if (ere !=0 ) result = bigRes;

  return result;
}

double IntegrandLambda_K(Aparams * app)
{
  
  oparams mopar = app->opar;
  double Llim = mopar.L;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  gsl_function F;

  F.function = &xIntegrandLambda_K;
  F.params = app;

  double result, error;
  gsl_set_error_handler_off();
  int ere = gsl_integration_qag (&F, -Llim, Llim, epsabs, epsilon, 1000, key, w,
				 &result,
				 &error);
  gsl_integration_workspace_free (w);

  if (ere !=0 ) result = bigRes;

  return result;
}

// Gamma and H
double IntegrandGH(Aparams * app)
{
  
  oparams mopar = app->opar;
  double Llim = mopar.L;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  gsl_function F;

  F.function = &xIntegrandGH;
  F.params = app;

  double result, error;
  gsl_set_error_handler_off();
  int ere = gsl_integration_qag (&F, -Llim, Llim, epsabs, epsilon, 1000, key, w,
				 &result,
				 &error);
  gsl_integration_workspace_free (w);

  if (ere !=0 ) result = bigRes;

  return result;
}

//======================================================================

double minim_Integrand(const gsl_vector *v, void *nullparams)
{
  std::cout << "Inside minim_Integration..." <<std::endl;
  struct oparams * mopar = static_cast<struct oparams *> (nullparams);
  double Llim = mopar->L;
  int ndim = mopar->ndim;

  double app[ndim];
  for (int ii=0;ii<ndim;ii++)
    {
      app[ii] = gsl_vector_get(v,ii);
    }

  Aparams  apar;
  apar.app = app;
  apar.opar = *mopar;

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  gsl_function F;

  F.function = &xIntegrand;
  F.params = &apar;

  double result, error;
  gsl_set_error_handler_off();
  int ere = gsl_integration_qag (&F, -Llim, Llim, epsabs, epsilon, 1000, key, w,
				 &result,
				 &error);
  if (ere !=0 ) result = bigRes;

  std::cout << " Coefficients" <<std::endl;
  for (int ii=0;ii<ndim;ii++)
    {
      std::cout << app[ii] << " ";
    }
  std::cout << std::endl;
  std::cout << "Results of integration " << result <<std::endl;
  gsl_integration_workspace_free (w);

  return result;
}

int main(int argc, char** argv)
{

#ifdef WITHGPERFTOOLS
  ProfilerStart("profile.log");
#endif

  clock_t timer; 
  timer = clock();

  std::string pname = argv[0];

  lineparser line;

  //  line.AddCommand('m',   "int"      ,   "1"  ,"model");
  //  line.AddCommand('f',   "string",   "table.dat", "CoefficientTable");
  line.AddCommand('c',   "string"   ,   "0110","EnergyTerms");
  line.AddCommand('L',   "double"   ,   "4.7","length"); 
  line.AddCommand('s',   "double"   ,   "1.e-2","Size for minimization to stop"); 
  line.AddCommand('o',   "string"   ,   "results.dat","output file");
  line.AddCommand('C',   "string"   ,   "coeffs.dat","output coefficient file");
  line.AddCommand('R',   "string"   ,   ""          ,"read coefficient file");
  line.AddCommand('e',   "double"   ,   "1.e-7","Relative epsilon"); 
  line.AddCommand('E',   "double"   ,   "1.e-3","Absolute epsilon"); 

  if (argc>1)
    {
      line.Interpret(argc,argv);
    }
  else
    {
      line.Print();
      error(argv[0]);
    }

  std::string cS = line.GetString('c');
  double Ls = line.GetDouble('L');

  // Set the absolute and relative errors for integrations.
  epsilon = line.GetDouble('e');
  epsabs = line.GetDouble('E');

  std::string read = line.GetString('R');
  
  std::cout << "Setting L to " << Ls << std::endl;
  if (cS.size() != 4)
    {
      std::cout << "Wrong EnergyTerms " << cS << std::endl;
      exit(1);
    }

  int c0,c2,c4,c6;
  std::cout << "Setting E0 term to ";
  if (cS[0] =='0')
    {
      c0 = 0;
      std::cout << " Off " << std::endl;
    }
  else if (cS[0]=='1')
    {
      c0 = 1;
      std::cout << " On " << std::endl;
    }
  else
    {
      std::cout << "Wrong EnergyTerms " << cS << std::endl;
      exit(1);
    }
  
  std::cout << "Setting E2 term to ";
  if (cS[1] =='0')
    {
      c2 = 0;
      std::cout << " Off " << std::endl;
    }
  else if (cS[1]=='1')
    {
      c2 = 1;
      std::cout << " On " << std::endl;
    }
  else
    {
      std::cout << "Wrong EnergyTerms " << cS << std::endl;
      exit(1);
    }
  
  std::cout << "Setting E4 term to ";
  if (cS[2] =='0')
    {
      c4 = 0;
      std::cout << " Off " << std::endl;
    }
  else if (cS[2]=='1')
    {
      c4 = 1;
      std::cout << " On " << std::endl;
    }
  else
    {
      std::cout << "Wrong EnergyTerms " << cS << std::endl;
      exit(1);
    }
  
  
  std::cout << "Setting E6 term to ";
  if (cS[3] =='0')
    {
      c6 = 0;
      std::cout << " Off " << std::endl;
    }
  else if (cS[3]=='1')
    {
      c6 = 1;
      std::cout << " On " << std::endl;
    }
  else
    {
      std::cout << "Wrong EnergyTerms " << cS << std::endl;
      exit(1);
    }
  


  // Order of the Fourier coefficients.
  int iorder = 6;

  // Number of Fourier coefficients as a function of the Order
  int nDimension[] = {2,8,16,31,32};
  // Index of the first beta coefficient 
  int nBeta[] ={1,4,9,18,24};
  
  // Number of fourier coefficients. 
  int ndim = nDimension[iorder-2];

  // Fits beta coefficient
  int iBeta = nBeta[iorder-2];

  std::cout << "#In Minim " << iorder << std::endl;
  std::cout << "# Number of dimensions " << ndim << std::endl;
  std::cout << "# First beta coeff. " << iBeta  << std::endl;
  if (iorder != 6)
    {
      std::cout << "Wrong iorder = "<< iorder << std::endl;
      exit(1);
    }

   double Llim = Ls;   // This redefines the limit of integration

   // Constants for the Model.
   struct oparams opar;
   opar.ndim = ndim; 
   opar.L = Llim;
   opar.c0 = c0; opar.c2 = c2; opar.c4 = c4; opar.c6 = c6; 


  double yparMin = 1.e-2; // Small value for the Fourier coeffs. 
  double yDelta = 1.e-1;  // initial delta step for the coeffs.

  // Value of size for which the convergence is assumed
  //  double sizeConvergence = 0.5e-2;
  double sizeConvergence = line.GetDouble('s');
  unsigned int Niter = 10000; // NMaximun number of iterations for convergence. 

  // Sine coefficients are started with -values
  // cosine coefficients are initially +
  double rset[ndim];
  for (int ii=0;ii< ndim; ii++)
    {
      rset[ii] = yparMin;
    }
  //#undef RANDOM 
#ifdef RANDOM
  std::cout << "I am randomizing the initial values "<< std::endl;
  const gsl_rng_type * uRan;
  gsl_rng * r;


  gsl_rng_env_setup();

  uRan = gsl_rng_default;
  r = gsl_rng_alloc (uRan);

  for (int i = 0; i < ndim; i++)
    {
      double u = gsl_rng_uniform (r);
      rset[i] = 2.0*(u-0.5)*rset[i];
      //printf ("%.5f\n", rset[i]);
    }

  gsl_rng_free (r);
#endif

  
  double ypar[ndim];
  for (int i = 0;i<ndim;i++)
    {
      ypar[i] = 0.0;
    }
  if (read == "" ) 
    {
    // Initial values of Fourier coefficients. All very small but the first
  // sin and cosine
      ypar[0] = 1.0; 
      ypar[iBeta] = -1.0;
    }
  else
  {
    if (iorder <= 2)
      {
	std::cout << " Can not read coefficients. " << std::endl;
	exit(1);
      }
    // If read is set, read the initial coefficients from previous run.
    ReadFile(read,ypar,ndim,iBeta,nDimension[iorder-3],nBeta[iorder-3]);
  }

  // Randomize the coeffs if needed.
  for (int i = 0;i<ndim;i++)
    {
      ypar[i] = ypar[i]+rset[i];
      //std::cout << i << " " << ypar[i] << std::endl; 
    }
  // Initial values of Fourier coefficients. All very small but the first
  // sin and cosine

  // Just for checking... 
  double x = 1.0;
  double y = 1.0;
  double z = 1.0;
  double E2p = E(ypar,x,y,z,opar);
  //std::cout <<" Value of E2 at r = (1,1,1)  "<< E2p << std::endl;
  
  struct Aparams apar; 
  apar.app = ypar;
  apar.opar = opar;

  double Etotal = Integrand(&apar);
  std::cout << " Quadrature using qag "<< std::endl;
  std::cout << "result  = " << Etotal << std::endl;

  std::cout << "Starting minimization..." << std::endl;

  const gsl_multimin_fminimizer_type *T =
    gsl_multimin_fminimizer_nmsimplex2;

  std::cout << "After minimizer..." << std::endl;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *xmia;
  gsl_multimin_function minex_func;

  std::cout << "After declarations..." << std::endl;

  size_t iter = 0;
  int status;
  double size;


  xmia = gsl_vector_alloc(ndim);
  ss = gsl_vector_alloc(ndim);
  for (int ii=0;ii<ndim; ii++)
    {
      gsl_vector_set(xmia,ii,ypar[ii]);
    }
  gsl_vector_set_all( ss, yDelta); //Initial step for the mimim grid

  std::cout<< "After initialization of parameters "  <<std::endl;

  minex_func.n = ndim;
  minex_func.f = minim_Integrand;
  minex_func.params = &opar;

  s = gsl_multimin_fminimizer_alloc (T, ndim);
  gsl_multimin_fminimizer_set (s, &minex_func, xmia, ss);

  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);

      if (status)
        break;

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, sizeConvergence);

      if (status == GSL_SUCCESS)
        {
          printf ("converged to minimum at\n");
        }

      printf ("Nint = %5d a0 = %10.3e a1 = %10.3e b0 = %10.3e  b1 = %10.3e  f() = %7.3f size = %.3f\n",
              int(iter),
              gsl_vector_get (s->x, 0),
              gsl_vector_get (s->x, 1),
              gsl_vector_get (s->x, iBeta),
              gsl_vector_get (s->x, iBeta+1),
              s->fval, size);
    }
  while (status == GSL_CONTINUE && iter < Niter);

  std::cout << "Final coefficients: " << std::endl;
  for (int ii=0;ii<ndim; ii++)
    {
      ypar[ii] = gsl_vector_get(s->x,ii);
      std::cout << ypar[ii] << "  ";
    }
  std::cout << std::endl;
  std::cout << "Final Enegy: " << s->fval << std::endl;
  double Efinal = s->fval;
  
  gsl_vector_free(xmia);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

  std::cout<< "Minimization completed. "  <<std::endl;
 
  struct Aparams finalpar; 
  finalpar.app = ypar;

  // Contributions to the energy
  opar.c0 = c0; opar.c2 = 0;opar.c4 = 0;opar.c6 = 0;
  finalpar.opar = opar;
  double E0 = Integrand(&finalpar);
  std::cout << "E0 = " << E0 << std::endl;

  opar.c0 = 0; opar.c2 = 1.0;opar.c4 = 0;opar.c6 = 0;
  finalpar.opar = opar;
  double E2 = Integrand(&finalpar);
  std::cout << "E2 = " << E2 << std::endl;

  opar.c0 = 0; opar.c2 = 0;opar.c4 = 1.0;opar.c6 = 0;
  finalpar.opar = opar;
  double E4 = Integrand(&finalpar);
  std::cout << "E4 = " << E4 << std::endl;

  opar.c0 = 0; opar.c2 = 0;opar.c4 = 0;opar.c6 = c6;
  finalpar.opar = opar;
  double E6 = Integrand(&finalpar);
  std::cout << "E6 = " << E6 << std::endl;

  // Contributions of the Kaon condensate potential
  opar.c0 = 1.0; opar.c2 = 0;opar.c4 = 0;opar.c6 = 0;
  finalpar.opar = opar;
  double A1 = IntegrandVk(&finalpar);
  std::cout << "A1 = " << A1 << std::endl;

  opar.c0 = 0; opar.c2 = 1.0;opar.c4 = 0;opar.c6 = 0;
  finalpar.opar = opar;
  double A2 = IntegrandVk(&finalpar);
  std::cout << "A2 = " << A2 << std::endl;

  opar.c0 = 0; opar.c2 = 0;opar.c4 = 1.0;opar.c6 = 0;
  finalpar.opar = opar;
  double A3 = IntegrandVk(&finalpar);
  std::cout << "A3 = " << A3 << std::endl;

  opar.c0 = 0; opar.c2 = 0;opar.c4 = 0;opar.c6 = c6;
  finalpar.opar = opar;
  double A4 = IntegrandVk(&finalpar);
  std::cout << "A4 = " << A4 << std::endl;
  
  // Contribution of the Isospin terms with Kaons
  opar.c0 = 1.0; opar.c2 = 0;opar.c4 = 0;opar.c6 = 0;
  finalpar.opar = opar;
  double L1 = IntegrandLambda_K(&finalpar);
  std::cout << "L1 = " << L1 << std::endl;

  opar.c0 = 0; opar.c2 = 1.0;opar.c4 = 0;opar.c6 = 0;
  finalpar.opar = opar;
  double L2 = IntegrandLambda_K(&finalpar);
  std::cout << "L2 = " << L2 << std::endl;

  opar.c0 = 0; opar.c2 = 0;opar.c4 = 1.0;opar.c6 = 0;
  finalpar.opar = opar;
  double L3 = IntegrandLambda_K(&finalpar);
  std::cout << "L3 = " << L3 << std::endl;
  
  // Gamma Contribution
  opar.c0 = 1.0; opar.c2 = 0;opar.c4 = 0;opar.c6 = 0;
  finalpar.opar = opar;
  double GammaL = IntegrandGH(&finalpar);
  std::cout << "Gamma = " << GammaL << std::endl;
  
  // H Contribution
  opar.c0 = 0; opar.c2 = 1.0;opar.c4 = 0;opar.c6 = 0;
  finalpar.opar = opar;
  double HL = IntegrandGH(&finalpar);
  std::cout << "H = " << HL << std::endl;
  
  // Topological number
  double Bn = IntegrandB(&finalpar);
  std::cout << "B = " << Bn << std::endl;

  // Isospin Contribution
  opar.c0 = c0; opar.c2 = c2; opar.c4 = c4; opar.c6 = c6;
  finalpar.opar = opar;
  double Isospin = IntegrandIso(&finalpar);
  std::cout << "Isospin = " << Isospin << std::endl;
  
  // Lambdaij Comprobation
  //double Lij = IntegrandLambda(&finalpar);
  //std::cout << "Lambda = " << Lij << std::endl;

  timer = clock() - timer;
  double seconds = (double(timer) )/CLOCKS_PER_SEC;

  std::string coeffile = line.GetString('C');
  std::ofstream cfile(coeffile.c_str());
  cfile << "# "<< Llim << " " << Efinal << " " << Bn << " " << Isospin << " " << A1 << " " << A2 << " " << A3 << " " << A4 << " " << L1 << " " << L2 << " " << L3 << " " << GammaL << " " << HL << std::endl;
	  
  for (int i=0;i<ndim;i++)
    {
      cfile << i << " " << ypar[i] << std::endl;
    }
  cfile.close();

  std::string ofile = line.GetString('o');
  std::ofstream file(ofile.c_str(),std::ios::app);
  file << Llim << " " << Efinal << " " << Bn << " " << Isospin << " " << A1 << " " << A2 << " " << A3 << " " << A4 << " " << L1 << " " << L2 << " " << L3 << " " << GammaL << " " << HL << std::endl;

  file.close();
  
  if (Llim == 4.7)
  {
  
   int Nx = 95;
   // Now write the Energy density for some points
   double xr[Nx] = {-4.7, -4.6, -4.5, -4.4, -4.3, -4.2, -4.1, -4. , -3.9, -3.8, -3.7,
       -3.6, -3.5, -3.4, -3.3, -3.2, -3.1, -3. , -2.9, -2.8, -2.7, -2.6,
       -2.5, -2.4, -2.3, -2.2, -2.1, -2. , -1.9, -1.8, -1.7, -1.6, -1.5,
       -1.4, -1.3, -1.2, -1.1, -1. , -0.9, -0.8, -0.7, -0.6, -0.5, -0.4,
       -0.3, -0.2, -0.1, -0. ,  0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,
        0.8,  0.9,  1. ,  1.1,  1.2,  1.3,  1.4,  1.5,  1.6,  1.7,  1.8,
        1.9,  2. ,  2.1,  2.2,  2.3,  2.4,  2.5,  2.6,  2.7,  2.8,  2.9,
        3. ,  3.1,  3.2,  3.3,  3.4,  3.5,  3.6,  3.7,  3.8,  3.9,  4. ,
        4.1,  4.2,  4.3,  4.4,  4.5,  4.6,  4.7};
       
   opar.c0 = c0; opar.c2 = c2; opar.c4 = c4; opar.c6 = c6;

   std::ofstream outfile1 ("ED_FCC_Half.dat");
   std::ofstream outfile2 ("pi1_FCC_Half.dat");
   std::ofstream outfile3 ("pi2_FCC_Half.dat");
   for (int i = 0; i<Nx; i++)
   {
    for (int j = 0; j<Nx; j++)
    {
      for (int k = 0; k<Nx; k++)
      {
       outfile1 << E(ypar, xr[i], xr[j], xr[k], opar) << "\t";
       outfile2 << B(ypar, xr[i], xr[j], xr[k],opar) << "\t";
       outfile3 << Iso(ypar, xr[i], xr[j], xr[k],opar) << "\t";
      }
    }
   }		

   outfile1.close();
   outfile2.close();
   outfile3.close();
  }

  /*
  if (Llim == 8.0)
  {
  
   opar.c0 = c0; opar.c2 = c2; opar.c4 = c4; opar.c6 = c6;
   int N = 41;
   
   double xg[N] = {-4.0, -3.9, -3.8, -3.7,
       -3.6, -3.5, -3.4, -3.3,
       -3.2, -3.1, -3.0, -2.9,
       -2.8, -2.7, -2.6, -2.5,
       -2.4, -2.3, -2.2, -2.1,
       -2.0, -1.9, -1.8, -1.7,
       -1.6, -1.5, -1.4, -1.3,
       -1.2, -1.1, -1.0, -0.9,
       -0.8, -0.7, -0.6, -0.5,
       -0.4, -0.3, -0.2, -0.1,
        0.0,  0.1,  0.2,  0.3,
        0.4,  0.5,  0.6,  0.7,
        0.8,  0.9,  1.0,  1.1,
        1.2,  1.3,  1.4,  1.5,
        1.6,  1.7,  1.8,  1.9,
        2.0,  2.1,  2.2,  2.3,
        2.4,  2.5,  2.6,  2.7,
        2.8,  2.9,  3.0,  3.1,
        3.2,  3.3,  3.4,  3.5,
        3.6,  3.7,  3.8,  3.9,
        4.0};
  
   std::ofstream outfileB ("rho.dat");
   
   for (int i = 0; i < N; i++)
   {
    for (int j = 0; j < N; j++)
    {
     for (int k = 0; k < N; k++)
     {
      outfileB << 0.5*B(ypar,xg[i],xg[j],xg[k],opar) - J03(ypar, xg[i], xg[j], xg[k], opar)/Isospin << "\t";
     }
    }
   }
   
   outfileB.close();
  }
  */
#ifdef WITHGPERFTOOLS
    ProfilerStop();
#endif
  return status;

}
