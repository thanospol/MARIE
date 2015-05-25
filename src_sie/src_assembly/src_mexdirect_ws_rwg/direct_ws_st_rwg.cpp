/**************************************************************************************************************************
           
		           DIRECT_WS_ST_RWG.cpp

Main body of the DIRECT EVALUATION method for the evaluation of the coincident 4-D
weakly singular integrals over planar triangular elements.

  Licensing: This code is distributed under the GNU LGPL license. 

  Modified:  21 October 2011

  Author:    Athanasios Polimeridis

  References

  A. G. Polimeridis and T. V. Yioultsis, “On the direct evaluation of weakly singular
  integrals in Galerkin mixed potential integral equation formulations,” IEEE Trans.
  Antennas Propag., vol. 56, no. 9, pp. 3011-3019, Sep. 2008.

  A. G. Polimeridis and J. R. Mosig, “Complete semi-analytical treatment of weakly
  singular integrals on planar triangles via the direct evaluation method,” Int. J.
  Numerical Methods Eng., vol. 83, pp. 1625-1650, 2010.

  A. G. Polimeridis, J. M. Tamayo, J. M. Rius and J. R. Mosig, “Fast and accurate
  computation of hyper-singular integrals in Galerkin surface integral equation
  formulations via the direct evaluation method,” IEEE Trans.
  Antennas Propag., vol. 59, no. 6, pp. 2329-2340, Jun. 2011.

  A. G. Polimeridis and J. R. Mosig, “On the direct evaluation of surface integral
  equation impedance matrix elements involving point singularities,” IEEE Antennas
  Wireless Propag. Lett., vol. 10, pp. 599-602, 2011.

  INPUT DATA
  r1,r2,r3 = point vectors of the triangular element's vertices
  Outer triangle P:(rp1,rp2,rp3)=(r1,r2,r3)
  Inner triangle Q:(rq1,rq2,rq3)=(r1,r2,r3)
  Np_1D = order of the Gauss-Legendre cubature for both dimensions of the remaining 1D smooth integral
  ko = wavenumber

  OUTPUT DATA
              I_DE(1)  = I_f1_f1
              I_DE(2)  = I_f1_f2
              I_DE(3)  = I_f1_f3
              I_DE(4)  = I_f2_f1
              I_DE(5)  = I_f2_f2
              I_DE(6)  = I_f2_f3
              I_DE(7)  = I_f3_f1
              I_DE(8)  = I_f3_f2
              I_DE(9)  = I_f3_f3

**************************************************************************************************************************/

#include "direct_ws_st_rwg.h"
#include <math.h>
using namespace std;

// ***********************************************************************
//			IMPLEMENTATION OF void DIRECT_ST_RWG
// ***********************************************************************

void direct_ws_st_rwg (const double r1[],const double r2[],const double r3[], const double ko, const int Np_1D, const double w[], const double z[], complex<double> I_DE[])
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

	complex<double> Isub[3][3][3];

	complex<double> I_simplex[10], I_vector[9], I_scalar;
	//

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

	 // for loops

	 for ( int column = 0 ; column <  3 ; column++ )
	 {
		 for ( int row = 0 ; row <  3 ; row++ )
		 {
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

					 PHI_a = PHI_functions(row, column, 1, PSI_1k, a_PSI_1k, ko);
					 PHI_e = PHI_functions(row, column, 5, PSI_1k, a_PSI_1k, ko);
					 PHI_f = PHI_functions(row, column, 6, PSI_1k, a_PSI_1k, ko);

					 F1 = PHI_a + PHI_e + PHI_f;       
					 Int1  = Int1 + w[kk] * F1;
					 // Int2,  pi/3 =< PSI <= 2pi/3
					 PSI_2a = M_PI / 3.0;
					 PSI_2b = 2.0 * M_PI / 3.0;
					 PSI_2k = ((PSI_2b - PSI_2a) * z[kk] + (PSI_2b + PSI_2a)) / 2.0;
					 a_PSI_2k = sqrt(acc * pow(cos(PSI_2k),2.0) - acs * cos(PSI_2k) * sin(PSI_2k) + ass * pow(sin(PSI_2k),2.0));

					 PHI_b = PHI_functions(row, column, 2, PSI_2k, a_PSI_2k, ko);
					 PHI_g = PHI_functions(row, column, 7, PSI_2k, a_PSI_2k, ko);

					 F2 = PHI_b + PHI_g ;       
					 Int2  = Int2 + w[kk] * F2;
					 // Int3,  2pi/3 =< PSI <= pi
					 PSI_3a = 2.0 * M_PI / 3.0;
					 PSI_3b = M_PI;
					 PSI_3k = ((PSI_3b - PSI_3a) * z[kk] + (PSI_3b + PSI_3a)) / 2.0;
					 a_PSI_3k = sqrt(acc * pow(cos(PSI_3k),2.0) - acs * cos(PSI_3k) * sin(PSI_3k) + ass * pow(sin(PSI_3k),2.0));

					 PHI_c = PHI_functions(row, column, 3, PSI_3k, a_PSI_3k, ko);
					 PHI_d = PHI_functions(row, column, 4, PSI_3k, a_PSI_3k, ko);
					 PHI_h = PHI_functions(row, column, 8, PSI_3k, a_PSI_3k, ko);

					 F3 = PHI_c + PHI_d + PHI_h ;       
					 Int3  = Int3 + w[kk] * F3;

				 }// for ( int kk = 0 ; kk <  Np_1D ; kk++ )

				 Int1 = ((PSI_1b - PSI_1a) / 2.0) * Int1;
				 Int2 = ((PSI_2b - PSI_2a) / 2.0) * Int2;
				 Int3 = ((PSI_3b - PSI_3a) / 2.0) * Int3;
				 //

				 Isub[row][column][m] = Int1 + Int2 + Int3;

			 }// for ( int m = 0 ; m <  3 ; m++ )

		 }// for ( int row = 0 ; row <  3 ; row++ )
	 }// for ( int column = 0 ; column <  3 ; column++ )

	// Post-processing

	DIRECT_post ( Isub, I_simplex );

	// 
	double L11  = vector_dot(l_1,l_1);
	double L22  = vector_dot(l_2,l_2);
	double L33  = vector_dot(l_3,l_3);
	//
	double L12  = vector_dot(l_1,l_2);
	double L13  = vector_dot(l_1,l_3);
	double L23  = vector_dot(l_2,l_3);
	//
	double L1   = sqrt(L11);
	double L2   = sqrt(L22);
	double L3   = sqrt(L33);
	//
	I_vector[0] = L22 * I_simplex[8] - 2 * L23 * I_simplex[5]     + L33 * I_simplex[4];
	I_vector[1] = L23 * I_simplex[6] - L12     * I_simplex[8]     - L33 * I_simplex[3] + L13 * I_simplex[5];
	I_vector[2] = L12 * I_simplex[7] - L22     * I_simplex[6]     - L13 * I_simplex[4] + L23 * I_simplex[3];
	I_vector[3] = L23 * I_simplex[2] - L33     * I_simplex[1]     - L12 * I_simplex[8] + L13 * I_simplex[7];
	I_vector[4] = L33 * I_simplex[0] - 2 * L13 * I_simplex[2]     + L11 * I_simplex[8];
	I_vector[5] = L13 * I_simplex[1] - L23     * I_simplex[0]     - L11 * I_simplex[7]  + L12 * I_simplex[6];
	I_vector[6] = L12 * I_simplex[5] - L13     * I_simplex[4]     - L22 * I_simplex[2]  + L23 * I_simplex[1];
	I_vector[7] = L13 * I_simplex[3] - L11     * I_simplex[5]     - L23 * I_simplex[0]  + L12 * I_simplex[2];
	I_vector[8] = L11 * I_simplex[4] - 2 * L12 * I_simplex[3]     + L22 * I_simplex[0];
	//
	I_scalar    = I_simplex[9];
	// Final output
	I_DE[0] = (L1*L1) / 12.0 * (Iunit * ko * I_vector[0] + 4.0 / (Iunit *ko) * I_scalar);
	I_DE[1] = (L1*L2) / 12.0 * (Iunit * ko * I_vector[1] + 4.0 / (Iunit *ko) * I_scalar);
	I_DE[2] = (L1*L3) / 12.0 * (Iunit * ko * I_vector[2] + 4.0 / (Iunit *ko) * I_scalar);
	I_DE[3] = (L2*L1) / 12.0 * (Iunit * ko * I_vector[3] + 4.0 / (Iunit *ko) * I_scalar);
	I_DE[4] = (L2*L2) / 12.0 * (Iunit * ko * I_vector[4] + 4.0 / (Iunit *ko) * I_scalar);
	I_DE[5] = (L2*L3) / 12.0 * (Iunit * ko * I_vector[5] + 4.0 / (Iunit *ko) * I_scalar);
	I_DE[6] = (L3*L1) / 12.0 * (Iunit * ko * I_vector[6] + 4.0 / (Iunit *ko) * I_scalar);
	I_DE[7] = (L3*L2) / 12.0 * (Iunit * ko * I_vector[7] + 4.0 / (Iunit *ko) * I_scalar);
	I_DE[8] = (L3*L3) / 12.0 * (Iunit * ko * I_vector[8] + 4.0 / (Iunit *ko) * I_scalar);
}

// ***********************************************************************
//			IMPLEMENTATION OF void DIRECT_post
// ***********************************************************************

void DIRECT_post ( complex<double> Isub[3][3][3], complex<double> I_simplex[] )
{
	complex<double> Ieta_xi[3][3];

	// Ieta_xi(1,1)
	Ieta_xi[0][0] = Isub[0][0][0] + Isub[0][0][1] + Isub[0][0][2];
	// Ieta_xi(1,2)
	complex<double> Ieta_xi_1_2_a =  Isub[0][1][0];
	complex<double> Ieta_xi_1_2_b =  0.5 * Isub[0][0][1] - 0.5 * Isub[0][1][1] - 0.5 * sqrt(3.0) * Isub[0][2][1];
	complex<double> Ieta_xi_1_2_c = -0.5 * Isub[0][0][2] - 0.5 * Isub[0][1][2] + 0.5 * sqrt(3.0) * Isub[0][2][2];

	Ieta_xi[0][1] = Ieta_xi_1_2_a + Ieta_xi_1_2_b + Ieta_xi_1_2_c;
	// Ieta_xi(1,3)
	complex<double> Ieta_xi_1_3_a =  Isub[0][2][0];
	complex<double> Ieta_xi_1_3_b =  0.5 * sqrt(3.0) * Isub[0][0][1] + 0.5 * sqrt(3.0) * Isub[0][1][1] - 0.5 * Isub[0][2][1];
	complex<double> Ieta_xi_1_3_c =  0.5 * sqrt(3.0) * Isub[0][0][2] - 0.5 * sqrt(3.0) * Isub[0][1][2] - 0.5 * Isub[0][2][2];

	Ieta_xi[0][2] = Ieta_xi_1_3_a + Ieta_xi_1_3_b + Ieta_xi_1_3_c;
	// Ieta_xi(2,1)
	complex<double> Ieta_xi_2_1_a =  Isub[1][0][0];
	complex<double> Ieta_xi_2_1_b =  0.5 * Isub[0][0][1] - 0.5 * Isub[1][0][1] - 0.5 * sqrt(3.0) * Isub[2][0][1];
	complex<double> Ieta_xi_2_1_c = -0.5 * Isub[0][0][2] - 0.5 * Isub[1][0][2] + 0.5 * sqrt(3.0) * Isub[2][0][2];
	
	Ieta_xi[1][0] = Ieta_xi_2_1_a + Ieta_xi_2_1_b + Ieta_xi_2_1_c;
	// Ieta_xi(2,2)
	complex<double> Ieta_xi_2_2_a =  Isub[1][1][0];
	complex<double> Ieta_xi_2_2_b =  0.25 *             Isub[0][0][1] - 0.25 *             Isub[0][1][1] - 0.25 * sqrt(3.0) * Isub[0][2][1]
	                               - 0.25 *             Isub[1][0][1] + 0.25 *             Isub[1][1][1] + 0.25 * sqrt(3.0) * Isub[1][2][1]
								   - 0.25 * sqrt(3.0) * Isub[2][0][1] + 0.25 * sqrt(3.0) * Isub[2][1][1] + 0.75 *             Isub[2][2][1];

	complex<double> Ieta_xi_2_2_c =  0.25 *             Isub[0][0][2] + 0.25 *             Isub[0][1][2] - 0.25 * sqrt(3.0) * Isub[0][2][2]
	                               + 0.25 *             Isub[1][0][2] + 0.25 *             Isub[1][1][2] - 0.25 * sqrt(3.0) * Isub[1][2][2]
								   - 0.25 * sqrt(3.0) * Isub[2][0][2] - 0.25 * sqrt(3.0) * Isub[2][1][2] + 0.75 *             Isub[2][2][2];

	Ieta_xi[1][1] = Ieta_xi_2_2_a + Ieta_xi_2_2_b + Ieta_xi_2_2_c;
	// Ieta_xi(2,3)
	complex<double> Ieta_xi_2_3_a =  Isub[1][2][0];
	complex<double> Ieta_xi_2_3_b =  0.25 * sqrt(3.0) * Isub[0][0][1] + 0.25 * sqrt(3.0) * Isub[0][1][1] - 0.25             * Isub[0][2][1]
	                               - 0.25 * sqrt(3.0) * Isub[1][0][1] - 0.25 * sqrt(3.0) * Isub[1][1][1] + 0.25 *             Isub[1][2][1]
								   - 0.75 *             Isub[2][0][1] - 0.75 *             Isub[2][1][1] + 0.25 * sqrt(3.0) * Isub[2][2][1];

	complex<double> Ieta_xi_2_3_c =- 0.25 * sqrt(3.0) * Isub[0][0][2] + 0.25 * sqrt(3.0) * Isub[0][1][2] + 0.25 *             Isub[0][2][2]
	                               - 0.25 * sqrt(3.0) * Isub[1][0][2] + 0.25 * sqrt(3.0) * Isub[1][1][2] + 0.25 *             Isub[1][2][2]
								   + 0.75 *             Isub[2][0][2] - 0.75 *             Isub[2][1][2] - 0.25 * sqrt(3.0) * Isub[2][2][2];

	Ieta_xi[1][2] = Ieta_xi_2_3_a + Ieta_xi_2_3_b + Ieta_xi_2_3_c;
	// Ieta_xi(3,1)
	complex<double> Ieta_xi_3_1_a =  Isub[2][0][0];
	complex<double> Ieta_xi_3_1_b =  0.5 * sqrt(3.0) * Isub[0][0][1] + 0.5 * sqrt(3.0) * Isub[1][0][1] - 0.5 * Isub[2][0][1];
	complex<double> Ieta_xi_3_1_c =  0.5 * sqrt(3.0) * Isub[0][0][2] - 0.5 * sqrt(3.0) * Isub[1][0][2] - 0.5 * Isub[2][0][2];
	
	Ieta_xi[2][0] = Ieta_xi_3_1_a + Ieta_xi_3_1_b + Ieta_xi_3_1_c;
	// Ieta_xi(3,2)
	complex<double> Ieta_xi_3_2_a =  Isub[2][1][0];
	complex<double> Ieta_xi_3_2_b =  0.25 * sqrt(3.0) * Isub[0][0][1] - 0.25 * sqrt(3.0) * Isub[0][1][1] - 0.75 *             Isub[0][2][1]
	                               + 0.25 * sqrt(3.0) * Isub[1][0][1] - 0.25 * sqrt(3.0) * Isub[1][1][1] - 0.75 *             Isub[1][2][1]
								   - 0.25 *             Isub[2][0][1] + 0.25 *             Isub[2][1][1] + 0.25 * sqrt(3.0) * Isub[2][2][1];

	complex<double> Ieta_xi_3_2_c =- 0.25 * sqrt(3.0) * Isub[0][0][2] - 0.25 * sqrt(3.0) * Isub[0][1][2] + 0.75 *             Isub[0][2][2]
	                               + 0.25 * sqrt(3.0) * Isub[1][0][2] + 0.25 * sqrt(3.0) * Isub[1][1][2] - 0.75 *             Isub[1][2][2]
								   + 0.25 *             Isub[2][0][2] + 0.25 *             Isub[2][1][2] - 0.25 * sqrt(3.0) * Isub[2][2][2];

	Ieta_xi[2][1] = Ieta_xi_3_2_a + Ieta_xi_3_2_b + Ieta_xi_3_2_c;
	// Ieta_xi(3,3)
	complex<double> Ieta_xi_3_3_a =  Isub[2][2][0];
	complex<double> Ieta_xi_3_3_b =  0.75 *             Isub[0][0][1] + 0.75 *             Isub[0][1][1] - 0.25 * sqrt(3.0) * Isub[0][2][1]
	                               + 0.75 *             Isub[1][0][1] + 0.75 *             Isub[1][1][1] - 0.25 * sqrt(3.0) * Isub[1][2][1]
								   - 0.25 * sqrt(3.0) * Isub[2][0][1] - 0.25 * sqrt(3.0) * Isub[2][1][1] + 0.25 *             Isub[2][2][1];

	complex<double> Ieta_xi_3_3_c =  0.75 *             Isub[0][0][2] - 0.75 *             Isub[0][1][2] - 0.25 * sqrt(3.0) * Isub[0][2][2]
	                               - 0.75 *             Isub[1][0][2] + 0.75 *             Isub[1][1][2] + 0.25 * sqrt(3.0) * Isub[1][2][2]
								   - 0.25 * sqrt(3.0) * Isub[2][0][2] + 0.25 * sqrt(3.0) * Isub[2][1][2] + 0.25 *             Isub[2][2][2];

	Ieta_xi[2][2] = Ieta_xi_3_3_a + Ieta_xi_3_3_b + Ieta_xi_3_3_c;
	// ***********************************************************************

	// 
	I_simplex[9] = Ieta_xi[0][0];

	// 
	complex<double> Ising_1_1 =  3.0 *       Ieta_xi[0][0] - 3.0 *       Ieta_xi[0][1] - sqrt(3.0) * Ieta_xi[0][2]
	                           - 3.0 *       Ieta_xi[1][0] + 3.0 *       Ieta_xi[1][1] + sqrt(3.0) * Ieta_xi[1][2]
							   - sqrt(3.0) * Ieta_xi[2][0] + sqrt(3.0) * Ieta_xi[2][1] +             Ieta_xi[2][2];

	I_simplex[0] = 1.0 / 12.0 * Ising_1_1;
	// 
	complex<double> Ising_1_2 =  3.0 *       Ieta_xi[0][0] + 3.0 *       Ieta_xi[0][1] - sqrt(3.0) * Ieta_xi[0][2]
	                           - 3.0 *       Ieta_xi[1][0] - 3.0 *       Ieta_xi[1][1] + sqrt(3.0) * Ieta_xi[1][2]
							   - sqrt(3.0) * Ieta_xi[2][0] - sqrt(3.0) * Ieta_xi[2][1] +             Ieta_xi[2][2];

	I_simplex[1] = 1.0 / 12.0 * Ising_1_2;
	// 
	complex<double> Ising_1_3 =  sqrt(3.0) * Ieta_xi[0][2] - sqrt(3.0) * Ieta_xi[1][2] -             Ieta_xi[2][2];

	I_simplex[2] = 1.0 / 6.0 * Ising_1_3;
	// 
	complex<double> Ising_2_1 =  3.0 *       Ieta_xi[0][0] - 3.0 *       Ieta_xi[0][1] - sqrt(3.0) * Ieta_xi[0][2]
	                           + 3.0 *       Ieta_xi[1][0] - 3.0 *       Ieta_xi[1][1] - sqrt(3.0) * Ieta_xi[1][2]
							   - sqrt(3.0) * Ieta_xi[2][0] + sqrt(3.0) * Ieta_xi[2][1] +             Ieta_xi[2][2];

	I_simplex[3] = 1.0 / 12.0 * Ising_2_1;
	// 
	complex<double> Ising_2_2 =  3.0 *       Ieta_xi[0][0] + 3.0 *       Ieta_xi[0][1] - sqrt(3.0) * Ieta_xi[0][2]
	                           + 3.0 *       Ieta_xi[1][0] + 3.0 *       Ieta_xi[1][1] - sqrt(3.0) * Ieta_xi[1][2]
							   - sqrt(3.0) * Ieta_xi[2][0] - sqrt(3.0) * Ieta_xi[2][1] +             Ieta_xi[2][2];

	I_simplex[4] = 1.0 / 12.0 * Ising_2_2;
	// 
	complex<double> Ising_2_3 =  sqrt(3.0) * Ieta_xi[0][2] + sqrt(3.0) * Ieta_xi[1][2] -             Ieta_xi[2][2];

	I_simplex[5] = 1.0 / 6.0 * Ising_2_3;
	// 
	complex<double> Ising_3_1 =  sqrt(3.0) * Ieta_xi[2][0] - sqrt(3.0) * Ieta_xi[2][1] -             Ieta_xi[2][2];

	I_simplex[6] = 1.0 / 6.0 * Ising_3_1;
	// 
	complex<double> Ising_3_2 =  sqrt(3.0) * Ieta_xi[2][0] + sqrt(3.0) * Ieta_xi[2][1] -             Ieta_xi[2][2];

	I_simplex[7] = 1.0 / 6.0 * Ising_3_2;
	// 
	I_simplex[8] = 1.0 / 3.0 * Ieta_xi[2][2];
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> PHI_functions
// ***********************************************************************

complex<double> PHI_functions ( int ROW, int COLUMN, int argument, double psi, double Apsi, const double ko )
{
	complex<double> PHI;
	complex<double> j   = Iunit;

	complex<double> PSI_1 ;
	complex<double> PSI_2 ;
	complex<double> PSI_3 ;
	complex<double> PSI_ ;

	complex<double> A = j * ko * Apsi;
	complex<double> B = cos(psi);
	complex<double> C = sin(psi) / sqrt(3.0);

	// Int_1_1 
	if (ROW == 0  && COLUMN == 0)
	{
		                                     
		PHI = PHI_1_1 (  A,  B,  C, argument,  psi,  Apsi,  ko );
		
	}// end if (ROW == 0  && COLUMN == 0)
	// Int_1_2
	else if (ROW == 0  && COLUMN == 1)
	{
		                                      
		PHI = PHI_1_2 (  A,  B,  C, argument,  psi,  Apsi,  ko );
		
	}// end if (ROW == 0  && COLUMN == 1)
	// Int_1_3
	else if (ROW == 0  && COLUMN == 2)
	{
		                                      
		PHI = PHI_1_3 (  A,  B,  C, argument,  psi,  Apsi,  ko );
		
	}// end if (ROW == 0  && COLUMN == 2)
	// Int_2_1
	else if (ROW == 1  && COLUMN == 0)
	{
		                                      
		PHI = PHI_2_1 (  A,  B,  C, argument,  psi,  Apsi,  ko );
		
	}// end if (ROW == 1  && COLUMN == 0)
	// Int_2_2
	else if (ROW == 1  && COLUMN == 1)
	{
		                                      
		PHI = PHI_2_2 (  A,  B,  C, argument,  psi,  Apsi,  ko );
		
	}// end if (ROW == 1  && COLUMN == 1)
	// Int_2_3
	else if (ROW == 1  && COLUMN == 2)
	{
		                                      
		PHI = PHI_2_3 (  A,  B,  C, argument,  psi,  Apsi,  ko );
		
	}// end if (ROW == 1  && COLUMN == 2)
	// Int_3_1
	else if (ROW == 2  && COLUMN == 0)
	{
		                                      
		PHI = PHI_3_1 (  A,  B,  C, argument,  psi,  Apsi,  ko );
		
	}// end if (ROW == 2  && COLUMN == 0)
	// Int_3_2
	else if (ROW == 2  && COLUMN == 1)
	{
		                                      
		PHI = PHI_3_2 (  A,  B,  C, argument,  psi,  Apsi,  ko );
		
	}// end if (ROW == 2  && COLUMN == 1)
	// Int_3_3
	else if (ROW == 2  && COLUMN == 2)
	{
		                                      
		PHI = PHI_3_3 (  A,  B,  C, argument,  psi,  Apsi,  ko );
		
	}// end if (ROW == 2  && COLUMN == 2)
	// Final Output
	return PHI;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> PHI_1_1
// ***********************************************************************

complex<double> PHI_1_1 ( complex<double> A, complex<double> B, complex<double> C, int argument, double psi, double Apsi, const double ko )
{
	complex<double> PHI_1_1;
	complex<double> j   = Iunit;

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

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> PHI_1_2
// ***********************************************************************

complex<double> PHI_1_2 ( complex<double> A, complex<double> B, complex<double> C, int argument, double psi, double Apsi, const double ko )
{
	complex<double> PHI_1_2;
	complex<double> j   = Iunit;

	complex<double> PSI_0 ;
	complex<double> PSI_1 ;
	complex<double> PSI_2_1 ;
	complex<double> PSI_2_2 ;
	complex<double> PSI_2_3 ;
	complex<double> PSI_2 ;
	complex<double> PSI_3 ;
	complex<double> PSI_4 ;
	complex<double> PSI_ ;

	double beta ;
	double gamma ;
	double epsilon ;
	double delta ;

		// Int_1_2                                      
		complex<double> D = sin(psi) / ( pow(j * ko,2.0) * pow(Apsi,2.0) );
		//
		switch(argument)
		{
			case 1:
				// PHI_a
				PSI_0 = 1.0 / (6.0 * B);
				PSI_1 = 1.0 / (2.0 * B);

				PSI_2_1 = 1.0;
				PSI_2_2 = -B / A * (1.0 - exp(-A / B) );
				PSI_2_3 = -1.0 / A * (B - (A + B) * exp(-A / B) );
				
				PSI_3  = 1.0 - B / A * (1.0 - exp(-A / B) );

				PSI_4  = 1.0 / 2.0 - B / A + (pow(B,2.0) / pow(A,2.0)) * (1.0 - exp(-A / B) );

				break;
			case 2:
				// PHI_b
				PSI_0 = 1.0 / (6.0 * C);
				PSI_1 = 1.0 / (2.0 * C);

				PSI_2_1 = 1.0;
				PSI_2_2 = -C / A * (1.0 - exp(-A / C) );
				PSI_2_3 = -1.0 / A * (C - (A + C) * exp(-A / C) );
				
				PSI_3  = 1.0 - C / A * (1.0 - exp(-A / C) );

				PSI_4  = 1.0 / 2.0 - C / A + (pow(C,2.0) / pow(A,2.0)) * (1.0 - exp(-A / C) );

				break;
			case 3:
				// PHI_c
				beta    = tan(M_PI - psi) / sqrt(3.0);
				gamma   = (1.0-beta) / (1.0+beta);

				PSI_0 = (1.0 / C) * ((1.0 / 2.0) * (1.0 - pow(gamma,2.0)) - (1.0 / 3.0) * (1.0 - pow(gamma,3.0) ) );
				PSI_1 = (1.0 / C) * (1.0 / 2.0 -(gamma - pow(gamma,2.0) / 2.0) );

				PSI_2_1 = 1.0 - gamma;
				PSI_2_2 = -C / A * (1.0 - exp(-A * (1.0 -gamma) / C) );
				PSI_2_3 = -1.0 /A * (C + (A * gamma - C - A) * exp(-A * (1.0 - gamma) / C) );
				
				PSI_3  = (1.0 - gamma) - C / A * (1.0 - exp(-A * (1.0-gamma) / C) );

				PSI_4  = 1.0 /2.0 * (1.0 - pow(gamma,2.0)) + C / pow(A,2.0) * (C - A + (A * gamma - C) * exp(-A * (1.0 - gamma) / C) );
								
				 break;
			case 4:
				// PHI_d
				beta    = tan(M_PI - psi) / sqrt(3.0);
				gamma   = (1.0-beta) / (1.0+beta);

				PSI_0 = -1.0 / B * (pow(gamma,2.0) / 2.0 + pow(gamma,3.0) / 3.0);
				PSI_1 = -1.0 / B * (gamma + pow(gamma,2.0) / 2.0 );

				PSI_2_1 = gamma;
				PSI_2_2 = -B / A * (exp(A * (1.0 + gamma) / B) - exp(A / B));
				PSI_2_3 = -1.0 / A * ((A - B) * exp(A / B) + (B - A - A * gamma) * exp(A * (1.0 + gamma) / B) );
				
				PSI_3  = gamma - B / A * (exp(A * (1.0 + gamma) / B) - exp(A /B) );

				PSI_4  = pow(gamma,2.0) / 2.0 - B / pow(A,2.0) * (B * exp(A / B) + (A * gamma - B) * exp(A * (1.0 + gamma) / B) );
								
				 break;
			case 5:
				// PHI_e
				epsilon = tan(psi) / sqrt(3.0);
				delta   = -(1-epsilon) / (1+epsilon);

				PSI_0 = 1.0 / B * (pow(delta,3.0) / 3.0 - pow(delta,2.0) / 2.0);
				PSI_1 = 1.0 / B * (pow(delta,2.0) / 2.0 - delta );

				PSI_2_1 = -delta;
				PSI_2_2 = -B / A * (exp(-A / B) - exp(-A * (1.0 - delta) / B) );
				PSI_2_3 = -1.0 / A * ((A + B) * exp(-A / B) + (A * delta - A - B) * exp(-A * (1.0 - delta) / B) );
				
				PSI_3  = -delta - B / A * (exp(-A / B) - exp(-A * (1.0 - delta) / B) );

				PSI_4  = -pow(delta,2.0) / 2.0 + B / pow(A,2.0) * (B * exp(-A / B) + (A * delta - B) * exp(-A * (1.0 - delta) / B) );
								
				 break;
			case 6:
				// PHI_f
				epsilon = tan(psi) / sqrt(3.0);
				delta   = -(1-epsilon) / (1+epsilon);

				PSI_0 = 1.0 / C * (1.0 / 2.0 * (pow(delta,2.0) - 1.0) + 1.0 /3.0 * (pow(delta,3.0) + 1.0) );
				PSI_1 = 1.0 / C * (delta + pow(delta,2.0) / 2.0 + 1.0 / 2.0);

				PSI_2_1 = 1.0 + delta;
				PSI_2_2 = -C / A * (1.0 - exp(-A * (1.0 + delta) / C) );
				PSI_2_3 = 1.0 / A * (( A + A * delta + C) * exp(-A * (1.0 + delta) / C) - C );
				
				PSI_3  = (1.0 + delta) - C / A * (1.0 - exp(-A * (1.0 + delta) / C) );

				PSI_4  = 1.0 / 2.0 * (pow(delta,2.0) - 1.0) + C / pow(A,2.0) * (A - C + (A * delta + C) * exp(-A * (1.0 + delta) / C) );
								
				 break;
			case 7:
				// PHI_g
				PSI_0 = -1.0 / (6.0 * C);
				PSI_1 = 1.0 / (2.0 * C);

				PSI_2_1 = 1.0 ;
				PSI_2_2 = -C / A * (1.0 - exp(-A / C) );
				PSI_2_3 = 1.0 / A * ((A + C) * exp(-A / C) - C);
				
				PSI_3  = 1.0 - C / A * (1.0 - exp(-A / C) );

				PSI_4  = -1.0 / 2.0 + C / pow(A,2.0) * (C * exp(-A / C) - C + A);
								
				 break;
			case 8:
				// PHI_h
				PSI_0 = 1.0 / (6.0 * B);
				PSI_1 = -1.0 / (2.0 * B);

				PSI_2_1 = 1.0 ;
				PSI_2_2 = -B / A * (exp(A / B) - 1.0);
				PSI_2_3 = 1.0 / A * ((A - B) * exp(A / B) + B);
				
				PSI_3  = 1.0 - B / A * (exp(A / B) - 1.0);

				PSI_4  = -1.0 / 2.0 - B / pow(A,2.0) * (A + B - B * exp(A / B) );
								
				 break;
		}// end switch argument

		PSI_2   = PSI_2_1 + PSI_2_2 + PSI_2_3;
		//
		PSI_ = j * ko *PSI_0 + (cos(psi) / Apsi) * PSI_1 - ((cos(psi) / (j * ko * pow(Apsi,2.0) ) ) * (PSI_2 + PSI_3) + (1.0 / Apsi) * PSI_4);
		//
		PHI_1_2    = D * PSI_;

	// Final Output
	return PHI_1_2;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> PHI_1_3
// ***********************************************************************

complex<double> PHI_1_3 ( complex<double> A, complex<double> B, complex<double> C, int argument, double psi, double Apsi, const double ko )
{
	complex<double> PHI_1_3;
	complex<double> j   = Iunit;

	complex<double> PSI_0 ;
	complex<double> PSI_1 ;
	complex<double> PSI_2 ;
	complex<double> PSI_ ;

	double beta ;
	double gamma ;
	double epsilon ;
	double delta ;

	// Int_1_3                                      
	complex<double> D = sin(psi) / ( pow(j * ko,2.0) * pow(Apsi,2.0) );
	//
	switch(argument)
		{
			case 1:
				// PHI_a
				PSI_0 = 1.0 / (3.0 * pow(B,2.0));
				PSI_1 = 1.0 / (2.0 * B);
				PSI_2 = 1.0 - B / A * (1.0 - exp(-A / B) );
				 break;
			case 2:
				// PHI_b
				PSI_0 =  1.0 / (3.0 * pow(C,2.0));
				PSI_1 = 1.0 / (2.0 * C);
				PSI_2 = 1.0 - C / A * (1.0 - exp(-A / C) );
				 break;
			case 3:
				// PHI_c
				beta    = tan(M_PI - psi) / sqrt(3.0);
				gamma   = (1.0-beta) / (1.0+beta);

				PSI_0 =  pow(1.0-gamma,3.0) / (3.0 * pow(C,2.0) );
				PSI_1 = (1.0 / C) * (1.0 / 2.0 - (gamma - pow(gamma,2.0) / 2.0 ) );
				PSI_2 = (1.0 - gamma) - C / A * (1.0 - exp(-A * (1.0 - gamma) / C) );
				 break;
			case 4:
				// PHI_d
				beta    = tan(M_PI - psi) / sqrt(3.0);
				gamma   = (1.0-beta) / (1.0+beta);

				PSI_0 = (pow(1.0 + gamma,3.0) - 1.0) / (3.0 * pow(B,2.0) );
				PSI_1 = -1.0 / B * (gamma + pow(gamma,2.0) / 2.0);
				PSI_2 = gamma - B / A * (exp(A * (1.0 + gamma) / B) - exp(A / B));
				 break;
			case 5:
				// PHI_e
				epsilon = tan(psi) / sqrt(3.0);
				delta   = -(1-epsilon) / (1+epsilon);

				PSI_0 =  (pow(1.0-delta,3.0) - 1.0) / (3.0 * pow(B,2.0) );
				PSI_1 = (1.0 / B) * ( pow(delta,2.0) / 2.0  - delta );
				PSI_2 = -delta - B / A * (exp(-A / B) - exp(-A * (1.0 - delta) / B));
				 break;
			case 6:
				// PHI_f
				epsilon = tan(psi) / sqrt(3.0);
				delta   = -(1-epsilon) / (1+epsilon);

				PSI_0 = pow(1.0 + delta,3.0) / (3.0 * pow(C,2.0));
				PSI_1 = 1.0 / C * (delta + pow(delta,2.0) / 2.0 + 1.0 / 2.0);
				PSI_2 = (1.0 + delta) - C / A * (1.0 - exp(-A * (1.0 + delta) / C) );
				 break;
			case 7:
				// PHI_g
				PSI_0 = 1.0 / (3.0 * pow(C,2.0) );
				PSI_1 = 1.0 / (2.0 * C);
				PSI_2 = 1.0 - C / A * (1.0 - exp(-A / C) );
				 break;
			case 8:
				// PHI_h
				PSI_0 = 1.0 / (3.0 * pow(B,2.0) );
				PSI_1 = -1.0 / (2.0 * B);
				PSI_2 = 1.0 - B / A * (exp(A / B) - 1.0);
				 break;
		}// end switch argument

		PSI_  = (j * ko * sin(psi) / 2.0) * PSI_0 - (sin(psi) / Apsi) * PSI_1 + (sin(psi) / (j * ko * pow(Apsi,2.0) ) ) * PSI_2;
		//
		PHI_1_3    = D * PSI_;

	// Final Output
	return PHI_1_3;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> PHI_2_1
// ***********************************************************************

complex<double> PHI_2_1 ( complex<double> A, complex<double> B, complex<double> C, int argument, double psi, double Apsi, const double ko )
{
	complex<double> PHI_2_1;
	complex<double> j   = Iunit;

	complex<double> PSI_1 ;
	complex<double> PSI_2 ;
	complex<double> PSI_3 ;

	double beta ;
	double gamma ;
	double epsilon ;
	double delta ;

	// Int_1_3                                      
	complex<double> D = sin(psi) / ( pow(j * ko,2.0) * pow(Apsi,3.0) );
	//
	switch(argument)
		{
			case 1:
				// PHI_a
				PSI_1 = A / (6.0 * B);
				PSI_2 = -1.0 / 2.0;
				PSI_3 = B / A - (pow(B,2.0) / pow(A,2.0) ) * (1.0 - exp(-A / B) );
				 break;
			case 2:
				// PHI_b
				PSI_1 = A / (6.0 * C);
				PSI_2 = -1.0 / 2.0;
				PSI_3 = C / A - (pow(C,2.0) / pow(A,2.0) ) * (1.0 - exp(-A / C) );
				 break;
			case 3:
				// PHI_c
				beta    = tan(M_PI - psi) / sqrt(3.0);
				gamma   = (1.0-beta) / (1.0+beta);

				PSI_1 = (A / C) * (1.0 / 2.0 * (1.0 - pow(gamma,2.0) ) - 1.0 / 3.0 * (1.0 - pow(gamma,3.0) ) );
				PSI_2 = -1.0 / 2.0 * (1.0 - pow(gamma,2.0) );
				PSI_3 = -(C / pow(A,2.0) ) * (C - A + (A * gamma - C) * exp(A * (gamma - 1.0) / C ) );
				 break;
			case 4:
				// PHI_d
				beta    = tan(M_PI - psi) / sqrt(3.0);
				gamma   = (1.0-beta) / (1.0+beta);

				PSI_1 = -(A / B) * (pow(gamma,2.0) / 2.0 + pow(gamma,3.0) / 3.0 );
				PSI_2 = -pow(gamma,2.0) / 2.0;
				PSI_3 = (B / pow(A,2.0)) * ((A * gamma - B) * exp(A * (1.0 + gamma) / B) + B * exp(A / B) );
				 break;
			case 5:
				// PHI_e
				epsilon = tan(psi) / sqrt(3.0);
				delta   = -(1-epsilon) / (1+epsilon);

				PSI_1 = A / B * (pow(delta,3.0) / 3.0 - pow(delta,2.0) / 2.0);
				PSI_2 = pow(delta,2.0) / 2.0;
				PSI_3 = -(B / pow(A,2.0)) * (B * exp(-A / B) + (A * delta - B) * exp(A * (delta - 1.0) / B ) );
				 break;
			case 6:
				// PHI_f
				epsilon = tan(psi) / sqrt(3.0);
				delta   = -(1-epsilon) / (1+epsilon);

				PSI_1 = (A / C ) * (1.0 / 2.0 *(pow(delta,2.0) - 1.0) + 1.0 / 3.0 * (pow(delta,3.0) + 1.0 ) );
				PSI_2 = -1.0 / 2.0 * (pow(delta,2.0) - 1.0);
				PSI_3 = -(C / pow(A,2.0) ) * (A - C + (A * delta + C) * exp(-A * (1.0 + delta) / C ) );
				 break;
			case 7:
				// PHI_g
				PSI_1 = -A / (6.0 * C);
				PSI_2 = 1.0 / 2.0;
				PSI_3 = -(C / pow(A,2.0)) * (C * exp(-A / C) - C + A);
				 break;
			case 8:
				// PHI_h
				PSI_1 = A / (6.0 * B);
				PSI_2 = 1.0 / 2.0;
				PSI_3 = (B / pow(A,2.0)) * (A + B - B * exp(A / B) );
				 break;
		}// end switch argument

		//
		PHI_2_1    = D * (PSI_1 + PSI_2 + PSI_3);

	// Final Output
	return PHI_2_1;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> PHI_2_2
// ***********************************************************************

complex<double> PHI_2_2 ( complex<double> A, complex<double> B, complex<double> C, int argument, double psi, double Apsi, const double ko )
{
	complex<double> PHI_2_2;
	complex<double> j   = Iunit;

	complex<double> PSI_0 ;
	complex<double> PSI_1 ;
	complex<double> PSI_2 ;
	complex<double> PSI_3 ;
	complex<double> PSI_3_1 ;
	complex<double> PSI_3_2 ;
	complex<double> PSI_4 ;
	complex<double> PSI_5 ;
	complex<double> PSI_ ;

	double beta ;
	double gamma ;
	double epsilon ;
	double delta ;

		// Int_2_2                                      
		complex<double> D = sin(psi) / ( pow(j * ko,2.0) * pow(Apsi,2.0) );
		//
		switch(argument)
		{
			case 1:
				// PHI_a
				PSI_0 = 1.0 / (6.0 * B);
				PSI_1 = 1.0 / (12.0 * B);
				PSI_2 = 1.0 / 2.0 - B / A + (pow(B,2.0) / pow(A,2.0)) * (1.0 - exp(-A / B) );

				PSI_5  = B / A - 2.0 * (pow(B,2.0) / pow(A,2.0)) + 2.0 * (pow(B,3.0) / pow(A,3.0)) * (1.0 - exp(-A / B) );
				
				PSI_3_1 = B / A - (pow(B,2.0) / pow(A,2.0)) * (1.0 - exp(-A / B));
				PSI_3_2 = PSI_5;
				PSI_3   = 1.0 / B * (PSI_3_1 - PSI_3_2);

				PSI_4  = 1.0 / 3.0 ;

				break;
			case 2:
				// PHI_b
				PSI_0 = 1.0 / (6.0 * C);
				PSI_1 = 1.0 / (12.0 * C);
				PSI_2 = 1.0 / 2.0 - C / A + (pow(C,2.0) / pow(A,2.0)) * (1.0 - exp(-A / C) );

				PSI_5  = C / A - 2.0 * (pow(C,2.0) / pow(A,2.0)) + 2.0 * (pow(C,3.0) / pow(A,3.0)) * (1.0 - exp(-A / C) );
				
				PSI_3_1 = C / A - (pow(C,2.0) / pow(A,2.0)) * (1.0 - exp(-A / C));
				PSI_3_2 = PSI_5;
				PSI_3   = 1.0 / C * (PSI_3_1 - PSI_3_2);

				PSI_4  = 1.0 / 3.0 ;

				break;
			case 3:
				// PHI_c
				beta    = tan(M_PI - psi) / sqrt(3.0);
				gamma   = (1.0-beta) / (1.0+beta);

				PSI_0 = (1.0 / C) * ((1.0 / 2.0) * (1.0 - pow(gamma,2.0)) - (1.0 / 3.0) * (1.0 - pow(gamma,3.0) ) );
				PSI_1 = (1.0 / C) * (1.0 / 12.0 - (pow(gamma,3.0) / 3.0 - pow(gamma,4.0) / 4.0) );
				PSI_2 = 1.0 / 2.0 * (1.0 - pow(gamma,2.0)) + C / pow(A,2.0) * (C - A + (A * gamma - C) * exp(-A * (1 - gamma) / C) );

				PSI_5  = C / A * (1.0 - pow(gamma,2.0) * exp(-A * (1.0 - gamma) / C) ) - 2.0 * (pow(C,2.0) / pow(A,2.0)) * (1.0- gamma * exp(-A * (1.0 - gamma) / C) ) + 2.0 * (pow(C,3.0) / pow(A,3.0)) * (1.0 - exp(-A * (1.0 - gamma) / C) );
				
				PSI_3_1 = -C / pow(A,2.0) * (C - A + (A * gamma - C) * exp(-A * (1.0 - gamma) / C) );
				PSI_3_2 = PSI_5;
				PSI_3   = 1.0 / C * (PSI_3_1 - PSI_3_2);

				PSI_4  = (1.0 - pow(gamma,3.0) ) / 3.0 ;
								
				 break;
			case 4:
				// PHI_d
				beta    = tan(M_PI - psi) / sqrt(3.0);
				gamma   = (1.0-beta) / (1.0+beta);

				PSI_0 = -1.0 / B * (pow(gamma,2.0) / 2.0 + pow(gamma,3.0) / 3.0);
				PSI_1 = -1.0 / B * (pow(gamma,3.0) / 3.0 + pow(gamma,4.0) / 4.0);
				PSI_2 = pow(gamma,2.0) / 2.0 - B / pow(A,2.0) * (B * exp(A / B) + (A * gamma - B) * exp(A * (1.0 + gamma) / B) );

				PSI_5  = B / A * pow(gamma,2.0) * exp(A * (1.0 + gamma) / B) - 2.0 * (pow(B,2.0) / pow(A,2.0)) * gamma * exp(A * (1.0 + gamma) / B) + 2.0 * (pow(B,3.0) / pow(A,3.0)) * (exp(A * (1.0 + gamma) / B) - exp(A / B) );
				
				PSI_3_1 = B / pow(A,2.0) * (B * exp(A / B) + (A * gamma - B) * exp(A * (1.0 + gamma) / B) );
				PSI_3_2 = PSI_5;
				PSI_3   = -1.0 / B * (PSI_3_1 + PSI_3_2);

				PSI_4  = pow(gamma,3.0) / 3.0;
								
				 break;
			case 5:
				// PHI_e
				epsilon = tan(psi) / sqrt(3.0);
				delta   = -(1-epsilon) / (1+epsilon);

				PSI_0 = 1.0 / B * (pow(delta,3.0) / 3.0 - pow(delta,2.0) / 2.0);
				PSI_1 = 1.0 / B * (pow(delta,4.0) / 4.0 - pow(delta,3.0) / 3.0);
				PSI_2 = -pow(delta,2.0) / 2.0 + B / pow(A,2.0) * (B * exp(-A / B) + (A * delta - B) * exp(-A * (1.0 - delta) / B) );

				PSI_5  = -B / A * pow(delta,2.0) * exp(-A * (1.0 - delta) / B) + 2.0 * (pow(B,2.0) / pow(A,2.0)) * delta * exp(-A * (1.0 - delta) / B) + 2.0 * (pow(B,3.0) / pow(A,3.0)) * (exp(-A / B) - exp(-A * (1.0 - delta) / B) );
				
				PSI_3_1 = -B / pow(A,2.0) * (B * exp(-A / B) + (A * delta - B) * exp(-A * (1.0 - delta) / B) );
				PSI_3_2 = PSI_5;
				PSI_3   = 1.0 / B * (PSI_3_1 - PSI_3_2);

				PSI_4  = -pow(delta,3.0) / 3.0;
								
				 break;
			case 6:
				// PHI_f
				epsilon = tan(psi) / sqrt(3.0);
				delta   = -(1-epsilon) / (1+epsilon);

				PSI_0 = 1.0 / C * (1.0 / 2.0 * (pow(delta,2.0) - 1.0) + 1.0 / 3.0 * (pow(delta,3.0) + 1.0) );
				PSI_1 = 1.0 / C * (pow(delta,3.0) / 3.0 + pow(delta,4.0) / 4.0 + 1.0 / 12.0);
				PSI_2 = 1.0 / 2.0 * (pow(delta,2.0) - 1.0) + C / pow(A,2.0) * (A - C + (A * delta + C) * exp(-A * (1.0 + delta) / C) );

				PSI_5  = C / A * (1.0 - pow(delta,2.0) * exp(-A * (1.0 + delta) / C)) - 2.0 * (pow(C,2.0) / pow(A,2.0)) * (delta * exp(-A * (1.0 + delta) / C) + 1.0) - 2.0 * (pow(C,3.0) / pow(A,3.0)) * (exp(-A * (1.0 + delta) / C) - 1.0);
								
				PSI_3_1 = -C / pow(A,2.0) * (A - C + (A * delta + C) * exp(-A * (1.0 + delta) / C) );
				PSI_3_2 = PSI_5;
				PSI_3   = 1.0 / C * (PSI_3_1 + PSI_3_2);

				PSI_4  = (pow(delta,3.0) + 1.0) / 3.0;
								
				 break;
			case 7:
				// PHI_g
				PSI_0 = -1.0 / (6.0 * C);
				PSI_1 = 1.0 / (12.0 * C);
				PSI_2 = -1.0 / 2.0 + C / pow(A,2.0) * (C * exp(-A / C) - C + A);

				PSI_5  = C / A - 2.0 * (pow(C,2.0) / pow(A,2.0)) - 2.0 * (pow(C,3.0) / pow(A,3.0)) * (exp(-A / C) - 1.0);
				
				PSI_3_1 = -C / pow(A,2.0) * (C * exp(-A / C) - C + A);
				PSI_3_2 = PSI_5;
				PSI_3   = 1.0 / C * (PSI_3_1 + PSI_3_2);

				PSI_4  = 1.0 / 3.0 ;
								
				 break;
			case 8:
				// PHI_h
				PSI_0 = 1.0 / (6.0 * B);
				PSI_1 = -1.0 / (12.0 * B);
				PSI_2 = -1.0 / 2.0 - B / pow(A,2.0) * ( A + B - B * exp(A / B) );

				PSI_5  = -B / A - 2.0 * pow(B,2.0) / pow(A,2.0) + 2.0 * (pow(B,3.0) / pow(A,3.0)) * (exp(A / B) - 1.0);
				
				PSI_3_1 = B / pow(A,2.0) * (A + B - B * exp(A / B) );
				PSI_3_2 = PSI_5;
				PSI_3   = -1.0 / B * (PSI_3_1 + PSI_3_2);

				PSI_4  = 1.0 / 3.0 ;
								
				 break;
		}// end switch argument
		
		PSI_ = (cos(psi) / Apsi) * PSI_0 + j * ko * PSI_1 - (2.0 * cos(psi) / (j * ko * pow(Apsi,2.0))) * PSI_2 + (cos(psi) / Apsi) * PSI_3 - (1 / Apsi) * (PSI_4 - PSI_5);
		//
		PHI_2_2    = D * PSI_;

	// Final Output
	return PHI_2_2;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> PHI_2_3
// ***********************************************************************

complex<double> PHI_2_3 ( complex<double> A, complex<double> B, complex<double> C, int argument, double psi, double Apsi, const double ko )
{
	complex<double> PHI_2_3;
	complex<double> j   = Iunit;

	complex<double> PSI_0 ;
	complex<double> PSI_1 ;
	complex<double> PSI_2 ;
	complex<double> PSI_ ;

	double beta ;
	double gamma ;
	double epsilon ;
	double delta ;

	// Int_1_3                                      
	complex<double> D = sin(psi) / ( pow(j * ko,2.0) * pow(Apsi,2.0) );
	//
	switch(argument)
		{
			case 1:
				// PHI_a
				PSI_0 = 1.0 / (12.0 * pow(B,2.0));
				PSI_1 = 1.0 / (6.0 * B);
				PSI_2 = 1.0/ 2.0 - B / A + (pow(B,2.0) / pow(A,2.0)) * (1.0 - exp(-A / B) );
				 break;
			case 2:
				// PHI_b
				PSI_0 =  1.0 / (12.0 * pow(C,2.0));
				PSI_1 = 1.0 / (6.0 * C);
				PSI_2 = 1.0 / 2.0 - C / A + (pow(C,2.0) / pow(A,2.0)) * (1.0 - exp(-A / C) );
				 break;
			case 3:
				// PHI_c
				beta    = tan(M_PI - psi) / sqrt(3.0);
				gamma   = (1.0-beta) / (1.0+beta);

				PSI_0 =  (1.0 / pow(C,2.0)) * (1.0 / 12.0 - ( pow(gamma,2.0) / 2.0 - 2.0 * pow(gamma,3.0) / 3.0 + pow(gamma,4.0) / 4.0) );
				PSI_1 = (1.0 / C) * ((1.0 / 2.0) * (1.0 - pow(gamma,2.0)) - (1.0 / 3.0) * (1.0 - pow(gamma,3.0) ) );
				PSI_2 = 1.0 / 2.0 * (1.0 - pow(gamma,2.0)) + C / pow(A,2.0) * (C - A + (A * gamma - C) * exp(-A * (1.0 - gamma) / C) );
				 break;
			case 4:
				// PHI_d
				beta    = tan(M_PI - psi) / sqrt(3.0);
				gamma   = (1.0-beta) / (1.0+beta);

				PSI_0 = (1.0 / pow(B,2.0)) * ( pow(gamma,2.0) / 2.0 + 2.0 * pow(gamma,3.0) / 3.0 + pow(gamma,4.0) / 4.0);
				PSI_1 = -1.0 / B * (pow(gamma,2.0) / 2.0 + pow(gamma,3.0) / 3.0);
				PSI_2 = pow(gamma,2.0) / 2.0 - B / pow(A,2.0) * (B * exp(A / B) + (A * gamma - B) * exp(A * (1.0 + gamma) / B) );
				 break;
			case 5:
				// PHI_e
				epsilon = tan(psi) / sqrt(3.0);
				delta   = -(1-epsilon) / (1+epsilon);

				PSI_0 =  -(1.0 / pow(B,2.0)) * (pow(delta,2.0) / 2.0 -2.0* pow(delta,3.0) / 3.0 + pow(delta,4.0) / 4.0);
				PSI_1 = 1.0 / B * (pow(delta,3.0) / 3.0 - pow(delta,2.0) / 2.0);
				PSI_2 = -pow(delta,2.0) / 2.0 + B / pow(A,2.0) * (B * exp(-A / B) + (A * delta - B) * exp(-A * (1.0 - delta) / B) );
				 break;
			case 6:
				// PHI_f
				epsilon = tan(psi) / sqrt(3.0);
				delta   = -(1-epsilon) / (1+epsilon);

				PSI_0 = (1.0 / pow(C,2.0)) * ((pow(delta,2.0) / 2.0 + 2.0* pow(delta,3.0) / 3.0 + pow(delta,4.0) / 4.0) - 1.0 / 12.0);
				PSI_1 = 1.0 / C * (1.0 / 2.0 * (pow(delta,2.0) - 1.0) + 1.0 / 3.0 * (pow(delta,3.0) + 1.0) );
				PSI_2 = 1.0 / 2.0 *(pow(delta,2.0) - 1.0) + C / pow(A,2.0) * (A - C + (A * delta + C) * exp(-A * (1.0 + delta) / C) );
				 break;
			case 7:
				// PHI_g
				PSI_0 = -1.0 / (12.0 * pow(C,2.0) );
				PSI_1 = -1.0 / (6.0 * C);
				PSI_2 = -1.0 / 2.0 + C / pow(A,2.0) * (C * exp(-A / C) - C + A);
				 break;
			case 8:
				// PHI_h
				PSI_0 = -1.0 / (12.0 * pow(B,2.0) );
				PSI_1 = 1.0 / (6.0 * B);
				PSI_2 = -1.0 / 2.0 - B / pow(A,2.0) * (A + B - B * exp(A / B) );
				 break;
		}// end switch argument

		PSI_  = (j * ko * sin(psi) / 2.0) * PSI_0 - (sin(psi) / Apsi) * PSI_1 + (sin(psi) / (j * ko * pow(Apsi,2.0)) ) * PSI_2;
		//
		PHI_2_3    = D * PSI_;

	// Final Output
	return PHI_2_3;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> PHI_3_1
// ***********************************************************************

complex<double> PHI_3_1 ( complex<double> A, complex<double> B, complex<double> C, int argument, double psi, double Apsi, const double ko )
{
	complex<double> PHI_3_1;
	complex<double> j   = Iunit;

	complex<double> PSI_2 ;
	complex<double> PSI_5_1 ;
	complex<double> PSI_5_2 ;
	complex<double> PSI_5_3 ;
	complex<double> PSI_ ;

	double beta ;
	double gamma ;
	double epsilon ;
	double delta ;

	// Int_1_3                                      
	complex<double> D = pow(sin(psi),2.0) / ( j * ko * pow(Apsi,2.0) );
	//
	switch(argument)
		{
			case 1:
				// PHI_a
				PSI_2 = 1.0 / (3.0 * pow(B,2.0));

				PSI_5_1 = 1.0 ;
				PSI_5_2 = -B / A * (1.0 - exp(-A / B) );
				PSI_5_3 = -1.0 / A * (B - (A + B) * exp(-A / B) );
				 break;
			case 2:
				// PHI_b
				PSI_2 = 1.0 / (3.0 * pow(C,2.0));

				PSI_5_1 = 1.0 ;
				PSI_5_2 = -C / A * (1.0 - exp(-A / C) );
				PSI_5_3 = -1.0 / A * (C - (A + C) * exp(-A / C) );
				 break;
			case 3:
				// PHI_c
				beta    = tan(M_PI - psi) / sqrt(3.0);
				gamma   = (1.0-beta) / (1.0+beta);

				PSI_2 = pow(1.0 - gamma,3.0) / (3.0 * pow(C,2.0));

				PSI_5_1 = 1.0 - gamma ;
				PSI_5_2 = -C / A * (1.0 - exp(-A * (1.0 - gamma) / C) );
				PSI_5_3 = -1.0 / A * (C + (A * gamma - C - A) * exp(-A * (1.0 - gamma) / C) );
				 break;
			case 4:
				// PHI_d
				beta    = tan(M_PI - psi) / sqrt(3.0);
				gamma   = (1.0-beta) / (1.0+beta);

				PSI_2 = (pow(1.0 + gamma,3.0) - 1.0 ) / (3.0 * pow(B,2.0) );

				PSI_5_1 = gamma ;
				PSI_5_2 = -B / A * (exp(A * (1.0 + gamma) / B) - exp(A / B) );
				PSI_5_3 = -1.0 / A * ((A - B) * exp(A / B) + (B - A - A * gamma) * exp(A * (1.0 + gamma) / B) );
				 break;
			case 5:
				// PHI_e
				epsilon = tan(psi) / sqrt(3.0);
				delta   = -(1-epsilon) / (1+epsilon);

				PSI_2 = (pow(1.0 - delta,3.0) - 1.0 ) / (3.0 * pow(B,2.0) );

				PSI_5_1 = -delta ;
				PSI_5_2 = -B / A * (exp(-A / B) - exp(-A * (1.0 - delta) / B) );
				PSI_5_3 = -1.0 / A * ((A + B) * exp(-A / B) + (A * delta - A - B) * exp(-A * (1.0 - delta) / B) );
				 break;
			case 6:
				// PHI_f
				epsilon = tan(psi) / sqrt(3.0);
				delta   = -(1-epsilon) / (1+epsilon);

				PSI_2 = pow(1.0 + delta,3.0) / (3.0 * pow(C,2.0) );

				PSI_5_1 = 1.0 + delta ;
				PSI_5_2 = -C / A * (1.0 - exp(-A * (1.0 + delta) / C) );
				PSI_5_3 = 1.0 / A * ((A + A * delta + C) * exp(-A * (1.0 + delta) / C) - C);
				 break;
			case 7:
				// PHI_g
				PSI_2 = 1.0 / (3.0 * pow(C,2.0) );

				PSI_5_1 = 1.0 ;
				PSI_5_2 = -C / A * (1.0 - exp(-A / C) );
				PSI_5_3 = 1.0 / A * ((A + C) * exp(-A / C) - C);
				 break;
			case 8:
				// PHI_h
				PSI_2 = 1.0 / (3.0 * pow(B,2.0) );

				PSI_5_1 = 1.0 ;
				PSI_5_2 = -B / A * (exp(A / B) - 1.0);
				PSI_5_3 = 1.0 / A * ((A - B) * exp(A / B) + B);
				 break;
		}// end switch argument

		PSI_  = PSI_2 / 2.0 - (PSI_5_1 + PSI_5_2 + PSI_5_3) / pow(j * ko * Apsi,2.0);
		//
		PHI_3_1    = D * PSI_;

	// Final Output
	return PHI_3_1;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> PHI_3_2
// ***********************************************************************

complex<double> PHI_3_2 ( complex<double> A, complex<double> B, complex<double> C, int argument, double psi, double Apsi, const double ko )
{
	complex<double> PHI_3_2;
	complex<double> j   = Iunit;

	complex<double> PSI_0 ;
	complex<double> PSI_1 ;
	complex<double> PSI_2_1 ;
	complex<double> PSI_2_2 ;
	complex<double> PSI_2_3 ;
	complex<double> PSI_2 ;
	complex<double> PSI_3 ;
	complex<double> PSI_3_1 ;
	complex<double> PSI_3_2 ;
	complex<double> PSI_3_3 ;
	complex<double> PSI_4 ;
	complex<double> PSI_5_1 ;
	complex<double> PSI_5_2 ;
	complex<double> PSI_5 ;
	complex<double> PSI_5_2_2 ;
	complex<double> PSI_ ;

	double beta ;
	double gamma ;
	double epsilon ;
	double delta ;

		// Int_1_2                                      
		complex<double> D = pow(sin(psi),2.0) / ( pow(j * ko,2.0) * pow(Apsi,2.0) );
		//
		switch(argument)
		{
			case 1:
				// PHI_a
				PSI_5_2_2 = B / A - 2.0 * pow(B,2.0) / pow(A,2.0) + 2.0 * pow(B,3.0) / pow(A,3.0) * (1.0 - exp(-A / B) );

				PSI_0 = 1.0 / (3.0 * pow(B,2.0));
				PSI_1 = 1.0 / (12.0 * pow(B,2.0));

				PSI_2_1 = 1.0;
				PSI_2_2 = -B / A * (1.0 - exp(-A / B) );
				PSI_2_3 = -1.0 / A * (B - (A + B) * exp(-A / B) );
				
				PSI_3_1  = B / A * (1.0 - exp(-A / B) );
				PSI_3_2  = B / A - pow(B,2.0) / pow(A,2.0) * (1.0 - exp(-A / B) );
				PSI_3_3  = PSI_5_2_2;

				PSI_3  = 1.0 / pow(B,2.0) * (PSI_3_1 - 2.0 * PSI_3_2 + PSI_3_3);

				PSI_4  = 1.0 / 2.0 - B / A + (pow(B,2.0) / pow(A,2.0)) * (1.0 - exp(-A / B) );

				PSI_5_1  = B / A - (pow(B,2.0) / pow(A,2.0)) * (1.0 - exp(-A / B) );
				PSI_5_2  = PSI_5_2_2;
				PSI_5  = 1.0 / B * (PSI_5_1 - PSI_5_2);

				break;
			case 2:
				// PHI_b
				PSI_5_2_2 = C / A - 2.0 * pow(C,2.0) / pow(A,2.0) + 2.0 * pow(C,3.0) / pow(A,3.0) * (1.0 - exp(-A / C) );

				PSI_0 = 1.0 / (3.0 * pow(C,2.0));
				PSI_1 = 1.0 / (12.0 * pow(C,2.0));

				PSI_2_1 = 1.0;
				PSI_2_2 = -C / A * (1.0 - exp(-A / C) );
				PSI_2_3 = -1.0 / A * (C - (A + C) * exp(-A / C) );
				
				PSI_3_1  = C / A * (1.0 - exp(-A / C) );
				PSI_3_2  = C / A - pow(C,2.0) / pow(A,2.0) * (1.0 - exp(-A / C) );
				PSI_3_3  = PSI_5_2_2;

				PSI_3  = 1.0 / pow(C,2.0) * (PSI_3_1 - 2.0 * PSI_3_2 + PSI_3_3);

				PSI_4  = 1.0 / 2.0 - C / A + (pow(C,2.0) / pow(A,2.0)) * (1.0 - exp(-A / C) );

				PSI_5_1  = C / A - (pow(C,2.0) / pow(A,2.0)) * (1.0 - exp(-A / C) );
				PSI_5_2  = PSI_5_2_2;
				PSI_5  = 1.0 / C * (PSI_5_1 - PSI_5_2);

				break;
			case 3:
				// PHI_c
				beta    = tan(M_PI - psi) / sqrt(3.0);
				gamma   = (1.0-beta) / (1.0+beta);

				PSI_5_2_2 = C / A * (1.0 - pow(gamma,2.0) * exp(-A * (1.0 - gamma) / C) ) - 2.0 * pow(C,2.0) / pow(A,2.0) * (1.0 - gamma * exp(-A * (1.0 - gamma) / C) ) + 2.0 * pow(C,3.0) / pow(A,3.0) * (1.0 - exp(-A * (1.0 - gamma) / C) );

				PSI_0 = pow(1.0 - gamma,3.0) / (3.0 * pow(C,2.0));
				PSI_1 = (1.0 / pow(C,2.0) ) * (1.0 / 12.0 - (pow(gamma,2.0) / 2.0 - 2.0 * pow(gamma,3.0) / 3.0 + pow(gamma,4.0) / 4.0) );

				PSI_2_1 = 1.0 - gamma;
				PSI_2_2 = -C / A * (1.0 - exp(-A * (1.0 - gamma) / C) );
				PSI_2_3 = -1.0 / A * (C + (A * gamma - C - A) * exp(-A * (1.0 - gamma) / C) );
				
				PSI_3_1  = C / A * (1.0 - exp(-A * (1.0 - gamma) / C) );
				PSI_3_2  = -C / pow(A,2.0) * (C - A + (A * gamma - C) * exp(-A * (1.0 - gamma) / C) );
				PSI_3_3  = PSI_5_2_2;

				PSI_3  = 1.0 / pow(C,2.0) * (PSI_3_1 - 2.0 * PSI_3_2 + PSI_3_3);

				PSI_4  = 1.0 / 2.0 * (1.0 - pow(gamma,2.0)) + C / pow(A,2.0) * (C - A + (A * gamma - C) * exp(-A * (1.0 - gamma) / C) );

				PSI_5_1  = -C / pow(A,2.0) * (C - A + (A * gamma - C) * exp(-A * (1.0 - gamma) / C) );
				PSI_5_2  = PSI_5_2_2;
				PSI_5  = 1.0 / C * (PSI_5_1 - PSI_5_2);
								
				 break;
			case 4:
				// PHI_d
				beta    = tan(M_PI - psi) / sqrt(3.0);
				gamma   = (1.0-beta) / (1.0+beta);

				PSI_5_2_2 = B / A * pow(gamma,2.0) * exp(A * (1.0 + gamma) / B) - 2.0 * pow(B,2.0) / pow(A,2.0) * gamma * exp(A * (1.0 + gamma) / B) + 2.0 * pow(B,3.0) / pow(A,3.0) * (exp(A * (1.0 + gamma) / B) - exp(A / B) );

				PSI_0 = (pow(1.0 + gamma,3.0) - 1.0) / (3.0 * pow(B,2.0));
				PSI_1 = (1.0 / pow(B,2.0) ) * ((pow(gamma,2.0) / 2.0 + 2.0 * pow(gamma,3.0) / 3.0 + pow(gamma,4.0) / 4.0) );

				PSI_2_1 = gamma;
				PSI_2_2 = -B / A * (exp(A * (1.0 + gamma) / B) - exp(A / B) );
				PSI_2_3 = -1.0 / A * ((A - B) * exp(A / B) + (B - A - A * gamma) * exp(A * (1.0 + gamma) / B) );
				
				PSI_3_1  = B / A * (exp(A * (1.0 + gamma) / B) - exp(A / B) );
				PSI_3_2  = B / pow(A,2.0) * (B * exp(A / B) + (A * gamma - B) * exp(A * (1.0 + gamma) / B) );
				PSI_3_3  = PSI_5_2_2;

				PSI_3  = 1.0 / pow(B,2.0) * (PSI_3_1 + 2.0 * PSI_3_2 + PSI_3_3);

				PSI_4  = pow(gamma,2.0) / 2.0 - B / pow(A,2.0) * (B * exp(A / B) + (A * gamma - B) * exp(A * (1.0 + gamma) / B) );

				PSI_5_1  =  B / pow(A,2.0) * (B * exp(A / B) + (A * gamma - B) * exp(A * (1.0 + gamma) / B) );
				PSI_5_2  = PSI_5_2_2;
				PSI_5  = -1.0 / B * (PSI_5_1 + PSI_5_2);
								
				 break;
			case 5:
				// PHI_e
				epsilon = tan(psi) / sqrt(3.0);
				delta   = -(1-epsilon) / (1+epsilon);

				PSI_5_2_2 = -B / A * pow(delta,2.0) * exp(-A * (1.0 - delta) / B) + 2.0 * pow(B,2.0) / pow(A,2.0) * delta * exp(-A * (1.0 - delta) / B) + 2.0 * pow(B,3.0) / pow(A,3.0) * (exp(-A / B) - exp(-A * (1.0 - delta) / B) );

				PSI_0 = (pow(1.0 - delta,3.0) - 1.0) / (3.0 * pow(B,2.0));
				PSI_1 = -(1.0 / pow(B,2.0) ) * ((pow(delta,2.0) / 2.0 - 2.0 * pow(delta,3.0) / 3.0 + pow(delta,4.0) / 4.0) );

				PSI_2_1 = -delta;
				PSI_2_2 = -B / A * (exp(- A / B) - exp(-A * (1.0 - delta) / B) );
				PSI_2_3 = -1.0 / A * ((A + B) * exp(-A / B) + (A * delta - A - B) * exp(-A * (1.0 - delta) / B) );
				
				PSI_3_1  = B / A * (exp(-A / B) - exp(-A * (1.0 - delta) / B) );
				PSI_3_2  = -B / pow(A,2.0) * (B * exp(-A / B) + (A * delta - B) * exp(-A * (1.0 - delta) / B) );
				PSI_3_3  = PSI_5_2_2;

				PSI_3  = 1.0 / pow(B,2.0) * (PSI_3_1 - 2.0 * PSI_3_2 + PSI_3_3);

				PSI_4  = -pow(delta,2.0) / 2.0 + B / pow(A,2.0) * (B * exp(-A / B) + (A * delta - B) * exp(-A * (1.0 - delta) / B) );

				PSI_5_1  = -B / pow(A,2.0) * (B * exp(-A / B) + (A * delta - B) * exp(-A * (1.0 - delta) / B) );
				PSI_5_2  = PSI_5_2_2;
				PSI_5  = 1.0 / B * (PSI_5_1 - PSI_5_2);
								
				 break;
			case 6:
				// PHI_f
				epsilon = tan(psi) / sqrt(3.0);
				delta   = -(1-epsilon) / (1+epsilon);

				PSI_5_2_2 = C / A * (1.0 - pow(delta,2.0) * exp(-A * (1.0 + delta) / C)) - 2.0 * pow(C,2.0) / pow(A,2.0) * (delta * exp(-A * (1.0 + delta) / C) + 1.0) - 2.0 * pow(C,3.0) / pow(A,3.0) * (exp(-A * (1.0 + delta) / C) - 1.0);

				PSI_0 = pow(1.0 + delta,3.0) / (3.0 * pow(C,2.0));
				PSI_1 = 1.0 / pow(C,2.0) * (pow(delta,2.0) / 2.0 + 2.0 * pow(delta,3.0) / 3.0 + pow(delta,4.0) / 4.0 - 1.0 / 12.0);

				PSI_2_1 = 1.0 + delta;
				PSI_2_2 = -C / A * (1.0 - exp(-A * (1.0 + delta) / C) );
				PSI_2_3 = 1.0 / A * ((A + A * delta + C) * exp(-A * (1.0 + delta) / C) - C);
				
				PSI_3_1  = C / A * (1.0 - exp(-A * (1.0 + delta) / C) );
				PSI_3_2  = -C / pow(A,2.0) * (A - C + (A * delta + C) * exp(-A * (1.0 + delta) / C) );
				PSI_3_3  = PSI_5_2_2;

				PSI_3  = 1.0 / pow(C,2.0) * (PSI_3_1 + 2.0 * PSI_3_2 + PSI_3_3);

				PSI_4  = 1.0 / 2.0 * (pow(delta,2.0) - 1.0) + C / pow(A,2.0) * (A - C + (A * delta + C) * exp(-A * (1.0 + delta) / C) );

				PSI_5_1  = -C / pow(A,2.0) * (A - C + (A * delta + C) * exp(-A * (1.0 + delta) / C) );
				PSI_5_2  = PSI_5_2_2;
				PSI_5  = 1.0 / C * (PSI_5_1 + PSI_5_2);
								
				 break;
			case 7:
				// PHI_g
				PSI_5_2_2 = C / A - 2.0 * pow(C,2.0) / pow(A,2.0) - 2.0 * pow(C,3.0) / pow(A,3.0) * (exp(-A / C) - 1.0);

				PSI_0 = 1.0 / (3.0 * pow(C,2.0));
				PSI_1 = -1.0 / (12.0 * pow(C,2.0));

				PSI_2_1 = 1.0 ;
				PSI_2_2 = -C / A * (1.0 - exp(-A / C) );
				PSI_2_3 = 1.0 / A * ((A + C) * exp(-A / C) - C);
				
				PSI_3_1  = C / A * (1.0 - exp(-A / C));
				PSI_3_2  = -C / pow(A,2.0) * (C * exp(-A / C) - C + A);
				PSI_3_3  = PSI_5_2_2;

				PSI_3  = 1.0 / pow(C,2.0) * (PSI_3_1 + 2.0 * PSI_3_2 + PSI_3_3);

				PSI_4  = -1.0 / 2.0 + C / pow(A,2.0) * (C * exp(-A / C) - C + A);

				PSI_5_1  = -C / pow(A,2.0) * (C * exp(-A / C) - C + A);
				PSI_5_2  = PSI_5_2_2;
				PSI_5  = 1.0 / C * (PSI_5_1 + PSI_5_2);
								
				 break;
			case 8:
				// PHI_h
				PSI_5_2_2 = -B / A - 2.0 * pow(B,2.0) / pow(A,2.0) + 2.0 * pow(B,3.0) / pow(A,3.0) * (exp(A / B) - 1.0);

				PSI_0 = 1.0 / (3.0 * pow(B,2.0));
				PSI_1 = -1.0 / (12.0 * pow(B,2.0));

				PSI_2_1 = 1.0 ;
				PSI_2_2 = -B / A * ( exp(A / B) - 1.0 );
				PSI_2_3 = 1.0 / A * ((A - B) * exp(A / B) + B);
				
				PSI_3_1  = B / A * (exp(A / B) -1.0);
				PSI_3_2  = B / pow(A,2.0) * (A +B -B * exp(A / B));
				PSI_3_3  = PSI_5_2_2;

				PSI_3  = 1.0 / pow(B,2.0) * (PSI_3_1 + 2.0 * PSI_3_2 + PSI_3_3);

				PSI_4  = -1.0 / 2.0 - B / pow(A,2.0) * (A + B - B * exp(A / B));

				PSI_5_1  = B / pow(A,2.0) * (A + B - B * exp(A / B));
				PSI_5_2  = PSI_5_2_2;
				PSI_5  = -1.0 / B * (PSI_5_1 + PSI_5_2);
								
				 break;
		}// end switch argument

		PSI_2 = PSI_2_1 + PSI_2_2 + PSI_2_3;
		//
		PSI_ = (cos(psi) / (2.0 * Apsi)) * PSI_0 + (j * ko / 2.0) * PSI_1 - (3.0 * cos(psi) / (pow(j * ko,2.0) * pow(Apsi,3.0)) ) * PSI_2 + (cos(psi) / Apsi) * PSI_3 - (1.0 / (j * ko * pow(Apsi,2.0)) ) * PSI_4 + (1.0 / Apsi) * PSI_5;
		//
		PHI_3_2    = D * PSI_;

	// Final Output
	return PHI_3_2;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> PHI_3_3
// ***********************************************************************

complex<double> PHI_3_3 ( complex<double> A, complex<double> B, complex<double> C, int argument, double psi, double Apsi, const double ko )
{
	complex<double> PHI_3_3;
	complex<double> j   = Iunit;

	complex<double> PSI_0 ;
	complex<double> PSI_1 ;
	complex<double> PSI_2_1 ;
	complex<double> PSI_2_2 ;
	complex<double> PSI_2_3 ;
	complex<double> PSI_2 ;
	complex<double> PSI_ ;

	double beta ;
	double gamma ;
	double epsilon ;
	double delta ;

	// Int_1_3                                      
	complex<double> D = pow(sin(psi),2.0) / ( pow(j * ko,2.0) * pow(Apsi,2.0) );
	//
	switch(argument)
		{
			case 1:
				// PHI_a
				PSI_0 = 1.0 / (4.0 * pow(B,3.0));
				PSI_1 = 1.0 / (3.0 * pow(B,2.0));

				PSI_2_1 = 1.0 ;
				PSI_2_2 = -B / A * (1.0 - exp(-A / B) );
				PSI_2_3 = -1.0 / A * (B - (A + B) * exp(-A / B) );
				 break;
			case 2:
				// PHI_b
				PSI_0 = 1.0 / (4.0 * pow(C,3.0));
				PSI_1 = 1.0 / (3.0 * pow(C,2.0));

				PSI_2_1 = 1.0 ;
				PSI_2_2 = -C / A * (1.0 - exp(-A / C) );
				PSI_2_3 = -1.0 / A * (C - (A + C) * exp(-A / C) );
				 break;
			case 3:
				// PHI_c
				beta    = tan(M_PI - psi) / sqrt(3.0);
				gamma   = (1.0-beta) / (1.0+beta);

				PSI_0 = pow(1.0 - gamma,4.0) / (4.0 * pow(C,3.0));
				PSI_1 = pow(1.0 - gamma,3.0) / (3.0 * pow(C,2.0));

				PSI_2_1 = 1.0 - gamma ;
				PSI_2_2 = -C / A * (1.0 - exp(-A * (1.0 - gamma) / C) );
				PSI_2_3 = -1.0 / A * (C + (A * gamma - C - A) * exp(-A * (1.0 - gamma) / C) );
				 break;
			case 4:
				// PHI_d
				beta    = tan(M_PI - psi) / sqrt(3.0);
				gamma   = (1.0-beta) / (1.0+beta);

				PSI_0 = (1.0 - pow(1.0 + gamma,4.0)) / (4.0 * pow(B,3.0));
				PSI_1 = ( pow(1.0 + gamma,3.0) - 1.0) / (3.0 * pow(B,2.0));

				PSI_2_1 = gamma ;
				PSI_2_2 = -B / A * (exp(A * (1.0 + gamma) / B) - exp(A / B) );
				PSI_2_3 = -1.0 / A * ((A - B) * exp(A / B) + (B - A - A * gamma) * exp(A * (1.0 + gamma) / B) );
				 break;
			case 5:
				// PHI_e
				epsilon = tan(psi) / sqrt(3.0);
				delta   = -(1-epsilon) / (1+epsilon);

				PSI_0 = (pow(1.0 - delta,4.0) - 1.0) / (4.0 * pow(B,3.0));
				PSI_1 = (pow(1.0 - delta,3.0) - 1.0) / (3.0 * pow(B,2.0));

				PSI_2_1 = -delta ;
				PSI_2_2 = -B / A * (exp(-A / B) - exp(-A * (1.0 - delta) / B) );
				PSI_2_3 = -1.0 / A * ((A + B) * exp(-A / B) + (A * delta - A - B) * exp(-A * (1.0 - delta) / B) );
				 break;
			case 6:
				// PHI_f
				epsilon = tan(psi) / sqrt(3.0);
				delta   = -(1-epsilon) / (1+epsilon);

				PSI_0 = pow(1.0 + delta,4.0) / (4.0 * pow(C,3.0));
				PSI_1 = pow(1.0 + delta,3.0) / (3.0 * pow(C,2.0));

				PSI_2_1 = 1.0 + delta ;
				PSI_2_2 = -C / A * (1.0 - exp(-A * (1.0 + delta) / C) );
				PSI_2_3 = 1.0 / A * ((A + A * delta + C) * exp(-A * (1.0 + delta) / C) - C);
				 break;
			case 7:
				// PHI_g
				PSI_0 =  1.0 / (4.0 * pow(C,3.0));
				PSI_1 =  1.0 / (3.0 * pow(C,2.0));

				PSI_2_1 = 1.0 ;
				PSI_2_2 = -C / A * (1.0 - exp(-A / C) );
				PSI_2_3 = 1.0 / A * ((A + C) * exp(-A / C) - C);
				 break;
			case 8:
				// PHI_h
				PSI_0 = - 1.0 / (4.0 * pow(B,3.0));
				PSI_1 =  1.0 / (3.0 * pow(B,2.0));

				PSI_2_1 = 1.0 ;
				PSI_2_2 = -B / A * (exp(A / B) - 1.0);
				PSI_2_3 = 1.0 / A * ((A - B) * exp(A / B) + B);
				 break;
		}// end switch argument

		PSI_2  = PSI_2_1+PSI_2_2+PSI_2_3;
		//
		PSI_ = (j * ko * sin(psi) / 3.0) * PSI_0 - (sin(psi) / (2.0 * Apsi)) *PSI_1 + (sin(psi) / (pow(j * ko,2.0) * pow(Apsi,3.0))) *PSI_2;
		//
		PHI_3_3    = D * PSI_;

	// Final Output
	return PHI_3_3;
}
