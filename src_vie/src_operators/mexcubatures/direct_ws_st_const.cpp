/**************************************************************************************************************************
           
		           DIRECT_WS_ST_RWG_const.cpp

Main body of the DIRECT EVALUATION method for the evaluation of the coincident 4-D
weakly singular integrals over planar triangular elements.

  Licensing: This code is distributed under the GNU LGPL license. 

  Modified:  14 March 2013

  Author:    Athanasios Polimeridis

  References

  A. G. Polimeridis and T. V. Yioultsis, �On the direct evaluation of weakly singular
  integrals in Galerkin mixed potential integral equation formulations,� IEEE Trans.
  Antennas Propag., vol. 56, no. 9, pp. 3011-3019, Sep. 2008.

  A. G. Polimeridis and J. R. Mosig, �Complete semi-analytical treatment of weakly
  singular integrals on planar triangles via the direct evaluation method,� Int. J.
  Numerical Methods Eng., vol. 83, pp. 1625-1650, 2010.

  A. G. Polimeridis, J. M. Tamayo, J. M. Rius and J. R. Mosig, �Fast and accurate
  computation of hyper-singular integrals in Galerkin surface integral equation
  formulations via the direct evaluation method,� IEEE Trans.
  Antennas Propag., vol. 59, no. 6, pp. 2329-2340, Jun. 2011.

  A. G. Polimeridis and J. R. Mosig, �On the direct evaluation of surface integral
  equation impedance matrix elements involving point singularities,� IEEE Antennas
  Wireless Propag. Lett., vol. 10, pp. 599-602, 2011.

  INPUT DATA
  r1,r2,r3 = point vectors of the triangular element's vertices
  Outer triangle P:(rp1,rp2,rp3)=(r1,r2,r3)
  Inner triangle Q:(rq1,rq2,rq3)=(r1,r2,r3)
  Np_1D = order of the Gauss-Legendre cubature for both dimensions of the remaining 1D smooth integral
  ko = wavenumber

  OUTPUT DATA
              I_DE  = I_const

**************************************************************************************************************************/

#include "direct_ws_st_const.h"
#include <math.h>
using namespace std;

// ***********************************************************************
//			IMPLEMENTATION OF void DIRECT_ST_RWG
// ***********************************************************************

complex<double> direct_ws_st_const (const double r1[],const double r2[],const double r3[], const double ko, const int Np_1D, const double z[], const double w[] )
{
	// ************************************************
	//			DECLARATION OF KEY VARIABLES
	// ************************************************

	// 1. Various

	double l_1[3], l_2[3], l_3[3], rp12[3], rp123[3];
	double Asing0, Asing1, Asing2, Asing3, Asing4, Asing5, Asing6, Asing7, Asing8;
	double Asing[3][3];
	double acc, acs, ass;

	complex<double> Int1, Int2, Int3;

	double PSI_1a, PSI_1b, PSI_1k, a_PSI_1k;
	double PSI_2a, PSI_2b, PSI_2k, a_PSI_2k;
	double PSI_3a, PSI_3b, PSI_3k, a_PSI_3k;

	complex<double> PHI_a, PHI_b, PHI_c, PHI_d, PHI_e, PHI_f, PHI_g, PHI_h, F1, F2, F3;

	complex<double> Isub[3];

	complex<double> I_simplex, I_scalar;
	// ************************************************
	//			MAIN CODE
	// ************************************************
	for (int i = 0; i < 3; i++)
	{
			rp12[i]    = r1[i] - r2[i];
			rp123[i]   = r1[i] + r2[i] - 2.0 * r3[i];
			//
			l_1[i]     = r2[i] - r3[i];
			l_2[i]     = r3[i] - r1[i];
			l_3[i]     = r1[i] - r2[i];
	}
    double Ap_[3];
    vector_cross(l_3, l_2, Ap_);
    double Ap = 1.0 / 2.0 * sqrt(vector_dot(Ap_, Ap_));
    double Jp = Ap / sqrt(3.0);
    double J = Jp * Jp;
	// Evaluate A parameters
	Asing0 = vector_dot(rp12,rp12) / 4.0;
	Asing1 = vector_dot(rp12,rp123)  / (2.0 * sqrt(3.0));
	Asing2 = vector_dot(rp123,rp123) / 12.0; 
	//
	Asing3 = (1.0/4.0)        *Asing0   -(sqrt(3.0)/4.0)  *Asing1      +(3.0/4.0)        *Asing2;
	Asing4 = (sqrt(3.0)/2.0)  *Asing0   -(1.0/2.0)        *Asing1      -(sqrt(3.0)/2.0)  *Asing2;
	Asing5 = (3.0/4.0)        *Asing0   +(sqrt(3.0)/4.0)  *Asing1      +(1.0/4.0)        *Asing2;
	//
	Asing6 = (1.0/4.0)        *Asing0   +(sqrt(3.0)/4.0)  *Asing1       +(3.0/4.0)       *Asing2;
	Asing7 = -(sqrt(3.0)/2.0) *Asing0   -(1.0/2.0)        *Asing1       +(sqrt(3.0)/2.0) *Asing2;
	Asing8 = (3.0/4.0)        *Asing0   -(sqrt(3.0)/4.0)  *Asing1       +(1.0/4.0)       *Asing2;
	//
	Asing[0][0] = Asing0;
	Asing[1][0] = Asing1;
	Asing[2][0] = Asing2;
	Asing[0][1] = Asing3;
	Asing[1][1] = Asing4;
	Asing[2][1] = Asing5;
	Asing[0][2] = Asing6;
	Asing[1][2] = Asing7;
	Asing[2][2] = Asing8;
	// Get the weights and abscissas for the 1-D quadrature
    for ( int m = 0 ; m <  3 ; m++ )
    {
        // Define the coefs of the appropriate subtriangle
        acc = Asing[0][m];
        acs = Asing[1][m];
        ass = Asing[2][m];
        // Gauss quadrature
        Int1 = 0.0;
        Int2 = 0.0;
        Int3 = 0.0;
        //
        for ( int kk = 0 ; kk <  Np_1D ; kk++ )
        {
            // Int1,  0 =< PSI <= pi/3
            PSI_1a = 0.0;
            PSI_1b = M_PI / 3.0;
            PSI_1k = ((PSI_1b - PSI_1a) * z[kk] + (PSI_1b + PSI_1a)) / 2.0;
            a_PSI_1k = sqrt(acc * pow(cos(PSI_1k),2.0) - acs * cos(PSI_1k) * sin(PSI_1k) + ass * pow(sin(PSI_1k),2.0));
            
            PHI_a = PHI_1_1( 1, PSI_1k, a_PSI_1k, ko);
            PHI_e = PHI_1_1( 5, PSI_1k, a_PSI_1k, ko);
            PHI_f = PHI_1_1( 6, PSI_1k, a_PSI_1k, ko);
            
            F1 = PHI_a + PHI_e + PHI_f;
            Int1  = Int1 + w[kk] * F1;
            // Int2,  pi/3 =< PSI <= 2pi/3
            PSI_2a = M_PI / 3.0;
            PSI_2b = 2.0 * M_PI / 3.0;
            PSI_2k = ((PSI_2b - PSI_2a) * z[kk] + (PSI_2b + PSI_2a)) / 2.0;
            a_PSI_2k = sqrt(acc * pow(cos(PSI_2k),2.0) - acs * cos(PSI_2k) * sin(PSI_2k) + ass * pow(sin(PSI_2k),2.0));
            
            PHI_b = PHI_1_1(2, PSI_2k, a_PSI_2k, ko);
            PHI_g = PHI_1_1(7, PSI_2k, a_PSI_2k, ko);
            
            F2 = PHI_b + PHI_g ;
            Int2  = Int2 + w[kk] * F2;
            // Int3,  2pi/3 =< PSI <= pi
            PSI_3a = 2.0 * M_PI / 3.0;
            PSI_3b = M_PI;
            PSI_3k = ((PSI_3b - PSI_3a) * z[kk] + (PSI_3b + PSI_3a)) / 2.0;
            a_PSI_3k = sqrt(acc * pow(cos(PSI_3k),2.0) - acs * cos(PSI_3k) * sin(PSI_3k) + ass * pow(sin(PSI_3k),2.0));
            
            PHI_c = PHI_1_1( 3, PSI_3k, a_PSI_3k, ko);
            PHI_d = PHI_1_1( 4, PSI_3k, a_PSI_3k, ko);
            PHI_h = PHI_1_1( 8, PSI_3k, a_PSI_3k, ko);
            
            F3 = PHI_c + PHI_d + PHI_h ;
            Int3  = Int3 + w[kk] * F3;
            
        }// for ( int kk = 0 ; kk <  Np_1D ; kk++ )
        
        Int1 = ((PSI_1b - PSI_1a) / 2.0) * Int1;
        Int2 = ((PSI_2b - PSI_2a) / 2.0) * Int2;
        Int3 = ((PSI_3b - PSI_3a) / 2.0) * Int3;
        //
        
        Isub[m] = Int1 + Int2 + Int3;
        
    }// for ( int m = 0 ; m <  3 ; m++ )
	// 
    I_simplex = Isub[0] + Isub[1] + Isub[2];
	//
	complex<double> I_DE    = J * I_simplex;
	// Final output
    return I_DE;
}

complex<double> PHI_1_1 ( int argument, double psi, double Apsi, const double ko )

//( complex<double> A, complex<double> B, complex<double> C, int argument, double psi, double Apsi, const double ko )

{

	complex<double> PHI_1_1;
	complex<double> j   = Iunit;
    
    complex<double> A = j * ko * Apsi;
	complex<double> B = cos(psi);
	complex<double> C = sin(psi) / sqrt(3.0);

	complex<double> PSI_1 ;
	complex<double> PSI_2 ;
	complex<double> PSI_3 ;
	complex<double> PSI_ ;

	double beta ;
	double gamma ;
	double epsilon ;
	double delta ;

    // Int_1_1
    complex<double> D = sin(psi) / ( pow(j*ko,2.0) * pow(Apsi,3.0) );
    //
    switch(argument)
    {
        case 1:
            // PHI_a
            PSI_1 = A / (2.0 * B);
            PSI_2 = -1.0;
            PSI_3 = (B / A) * (1.0 - exp(-A/B));
            break;
        case 2:
            // PHI_b
            PSI_1 =  A / (2.0 * C);
            PSI_2 = -1.0;
            PSI_3 = (C / A) * (1.0 - exp(-A/C));
            break;
        case 3:
            // PHI_c
            beta    = tan(M_PI - psi) / sqrt(3.0);
            gamma   = (1.0-beta) / (1.0+beta);
            
            PSI_1 =  (A/C) * (1.0/2.0 -(gamma-pow(gamma,2.0) / 2.0) );
            PSI_2 = - (1.0 - gamma);
            PSI_3 = (C / A) * (1.0 - exp(-A * (1-gamma) / C) );
            break;
        case 4:
            // PHI_d
            beta    = tan(M_PI - psi) / sqrt(3.0);
            gamma   = (1.0-beta) / (1.0+beta);
            
            PSI_1 = -A * (gamma + pow(gamma,2.0) / 2.0) / B;
            PSI_2 = -gamma;
            PSI_3 = (B / A) * (exp(A * (1+gamma) / B) - exp(A/B));
            break;
        case 5:
            // PHI_e
            epsilon = tan(psi) / sqrt(3.0);
            delta   = -(1-epsilon) / (1+epsilon);
            
            PSI_1 = A * (pow(delta,2.0) / 2.0 - delta) / B;
            PSI_2 = delta;
            PSI_3 = (B / A) * (exp(-A/B) - exp(-A * (1-delta) / B) );
            break;
        case 6:
            // PHI_f
            epsilon = tan(psi) / sqrt(3.0);
            delta   = -(1-epsilon) / (1+epsilon);
            
            PSI_1 = (A / C) * (delta + pow(delta,2.0) / 2.0 + 1.0/2.0);
            PSI_2 = -(1.0 + delta);
            PSI_3 = (C / A) * (1.0 - exp(-A * (1 + delta) / C));
            break;
        case 7:
            // PHI_g
            PSI_1 = A / (2.0 * C);
            PSI_2 = -1.0;
            PSI_3 = (C / A) * (1.0 - exp(-A / C) );
            break;
        case 8:
            // PHI_h
            PSI_1 = -A / (2.0 * B);
            PSI_2 = -1.0;
            PSI_3 = (B / A) * (exp(A/B) - 1.0);
            break;
    }// end switch argument
    
    PSI_  = PSI_1 + PSI_2 + PSI_3;
    //
    PHI_1_1    = D * PSI_;

	// Final Output
	return PHI_1_1;
}