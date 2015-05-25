/**************************************************************************************************************************
           
		           DIRECT_WS_VA_RWG.cpp

Main body of the DIRECT EVALUATION method for the evaluation of the edge adjacent 4-D
weakly singular integrals over planar triangular elements.

  Licensing: This code is distributed under the GNU LGPL license. 

  Modified:  24 October 2011

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
  r1,r2,r3,r4,r5 = point vectors of the triangular element's vertices
  Outer triangle P:(rp1,rp2,rp3)=(r1,r2,r3)
  Inner triangle Q:(rq1,rq2,rq3)=(r1,r2,r5)
  N_theta_p = order of the GL quadrature for the 1-D smooth integral over theta_p
  N_theta_q = order of the GL quadrature for the 1-D smooth integral over theta_q
  N_psi     = order of the GL quadrature for the 1-D smooth integral over Psi
  ko = wavenumber

  OUTPUT DATA
              I_DE[0]  = I_f1_f1
              I_DE[1]  = I_f1_f2
              I_DE[2]  = I_f1_f3
              I_DE[3]  = I_f2_f1
              I_DE[4]  = I_f2_f2
              I_DE[5]  = I_f2_f3
              I_DE[6]  = I_f3_f1
              I_DE[7]  = I_f3_f2
              I_DE[8]  = I_f3_f3

**************************************************************************************************************************/

#include "direct_ws_va_rwg.h"

// ***********************************************************************
//			IMPLEMENTATION OF void DIRECT
// ***********************************************************************

void direct_ws_va_rwg (const double r1[],const double r2[],const double r3[],const double r4[], const double r5[], const double ko, const int N_theta_p, const int N_theta_q, const int N_psi, const double w_theta_p[], const double z_theta_p[], const double w_theta_q[], const double z_theta_q[], const double w_psi[], const double z_psi[], complex<double> I_DE[] )
{
	// ************************************************
	//			DECLARATION OF KEY VARIABLES
	// ************************************************

	// 1. Various

	double r21[3], r31[3], r32[3], r41[3], r54[3], r51[3];
	//
	complex<double> Ising_const, Ising_f1_f1, Ising_f1_f2, Ising_f1_f3, Ising_f2_f1, Ising_f2_f2, Ising_f2_f3, Ising_f3_f1, Ising_f3_f2, Ising_f3_f3;
	//
	double theta_p_A, theta_p_B, THETA_p, J_theta_p, Lp;
	double theta_q_A, theta_q_B, THETA_q, J_theta_q, Lq;

	double b_v0[3], b_v1[3];

	double alpha_vo[3], alpha_v1[3], alpha_v2[3], alpha_v3[3];

	double a_v;
	double b_v;
	double c_v;

	
	double psi_1_A, psi_1_B, PSI_1, J_psi_1, psi_2_A, psi_2_B, PSI_2, J_psi_2;

	
	double B_1;
	double B_2;

	double L1, L2;

	complex<double> K_1[4], K_2[4];

	complex<double> Omega_f1_f1_1, Omega_f1_f2_1, Omega_f1_f3_1, Omega_f2_f1_1, Omega_f2_f2_1, Omega_f2_f3_1, Omega_f3_f1_1, Omega_f3_f2_1, Omega_f3_f3_1, Omega_const_1;
	complex<double> Omega_f1_f1_2, Omega_f1_f2_2, Omega_f1_f3_2, Omega_f2_f1_2, Omega_f2_f2_2, Omega_f2_f3_2, Omega_f3_f1_2, Omega_f3_f2_2, Omega_f3_f3_2, Omega_const_2;
	//
	double WPSI, WTHETA_p, WTHETA_q;
	//
	complex<double> I_theta_p_1_1, I_theta_p_1_2, I_theta_p_1_3, I_theta_p_2_1, I_theta_p_2_2, I_theta_p_2_3, I_theta_p_3_1, I_theta_p_3_2, I_theta_p_3_3, I_theta_p_const; 
	complex<double> I_theta_q_1_1, I_theta_q_1_2, I_theta_q_1_3, I_theta_q_2_1, I_theta_q_2_2, I_theta_q_2_3, I_theta_q_3_1, I_theta_q_3_2, I_theta_q_3_3, I_theta_q_const; 
	complex<double> I_psi_1_1, I_psi_1_2, I_psi_1_3, I_psi_2_1, I_psi_2_2, I_psi_2_3, I_psi_3_1, I_psi_3_2, I_psi_3_3, I_psi_const; 
	//

	// 2. Coefficients' parameters

	complex<double> coef_const[1];
	complex<double> coef_f1_f1[4];
	complex<double> coef_f1_f2[6];
	complex<double> coef_f1_f3[6];
	complex<double> coef_f2_f1[6];
	complex<double> coef_f2_f2[9];
	complex<double> coef_f2_f3[9];
	complex<double> coef_f3_f1[6];
	complex<double> coef_f3_f2[9];
	complex<double> coef_f3_f3[9];

	// ************************************************
	//			MAIN CODE
	// ************************************************
	// Evaluate alpha, beta, gamma r21, r31, r32, r41, r42 parameters
	for (int i = 0; i < 3; i++)
	{
			alpha_vo[i] = (r2[i] - r1[i]) / 2;
			alpha_v1[i] = (2 * r3[i] - r1[i] - r2[i]) / (2.0 * sqrt(3.0));
			alpha_v2[i] = (r4[i] - r1[i]) / 2;
			alpha_v3[i] = (2 * r5[i] - r1[i] - r4[i]) / (2.0 * sqrt(3.0));
			r21[i]   = r2[i] - r1[i];
			r31[i]   = r3[i] - r1[i];
			r32[i]   = r3[i] - r2[i];
			r41[i]   = r4[i] - r1[i];
	        r51[i]   = r5[i] - r1[i];
			r54[i]   = r5[i] - r4[i];
	}
	// Get the coefficients
	coefficients_const ( r1,  r2,  r3,  r4, r5,  ko,  coef_const );
	coefficients_f1_f1 ( r1,  r2,  r3,  r4, r5,  ko,  coef_f1_f1 );
	coefficients_f1_f2 ( r1,  r2,  r3,  r4, r5,  ko,  coef_f1_f2 );
	coefficients_f1_f3 ( r1,  r2,  r3,  r4, r5,  ko,  coef_f1_f3 );
	coefficients_f2_f1 ( r1,  r2,  r3,  r4, r5,  ko,  coef_f2_f1 );
	coefficients_f2_f2 ( r1,  r2,  r3,  r4, r5,  ko,  coef_f2_f2 );
	coefficients_f2_f3 ( r1,  r2,  r3,  r4, r5,  ko,  coef_f2_f3 );
	coefficients_f3_f1 ( r1,  r2,  r3,  r4, r5,  ko,  coef_f3_f1 );
	coefficients_f3_f2 ( r1,  r2,  r3,  r4, r5,  ko,  coef_f3_f2 );
	coefficients_f3_f3 ( r1,  r2,  r3,  r4, r5,  ko,  coef_f3_f3 );
     // Initialization of I_theta_p
		I_theta_p_1_1 = 0.0;
		I_theta_p_1_2 = 0.0;
		I_theta_p_1_3 = 0.0;
		I_theta_p_2_1 = 0.0;
		I_theta_p_2_2 = 0.0;
		I_theta_p_2_3 = 0.0;
		I_theta_p_3_1 = 0.0;
		I_theta_p_3_2 = 0.0;
		I_theta_p_3_3 = 0.0;
		I_theta_p_const = 0.0;
	 //
	 for ( int n_theta_p = 0 ; n_theta_p <  N_theta_p ; n_theta_p++ )
	 {
		theta_p_A = 0.0;
        theta_p_B = M_PI / 3.0;
        //
        THETA_p = ((theta_p_B - theta_p_A) / 2.0) * z_theta_p[n_theta_p] + (theta_p_B + theta_p_A) / 2.0;
        //
        J_theta_p = (theta_p_B - theta_p_A) / 2.0;
        //
        Lp = (2.0 * sqrt(3.0) ) / (sin(THETA_p) + sqrt(3.0) * cos(THETA_p) );
        //
			for (int i = 0; i < 3; i++)
				{
					b_v0[i]   = alpha_vo[i] * cos(THETA_p) + alpha_v1[i] * sin(THETA_p);
				}
		// Initialization of I_theta_q
		I_theta_q_1_1 = 0.0;
		I_theta_q_1_2 = 0.0;
		I_theta_q_1_3 = 0.0;
		I_theta_q_2_1 = 0.0;
		I_theta_q_2_2 = 0.0;
		I_theta_q_2_3 = 0.0;
		I_theta_q_3_1 = 0.0;
		I_theta_q_3_2 = 0.0;
		I_theta_q_3_3 = 0.0;
		I_theta_q_const = 0.0;
		//
		 for ( int n_theta_q = 0 ; n_theta_q <  N_theta_q ; n_theta_q++ )
		 {
			theta_q_A = 0.0;
			theta_q_B = M_PI / 3.0;
			//
			THETA_q = ((theta_q_B - theta_q_A) / 2.0) * z_theta_q[n_theta_q] + (theta_q_B + theta_q_A) / 2.0;
			//
			J_theta_q = (theta_q_B - theta_q_A) / 2.0;
			//
			Lq = (2.0 * sqrt(3.0) ) / (sin(THETA_q) + sqrt(3.0) * cos(THETA_q) );
			//
			for (int i = 0; i < 3; i++)
				{
					b_v1[i]   = alpha_v2[i] * cos(THETA_q) + alpha_v3[i] * sin(THETA_q);
				}
			//
			a_v = vector_dot(b_v0,b_v0);
			b_v = -2.0 * vector_dot(b_v0,b_v1);
			c_v = vector_dot(b_v1,b_v1);
			// Initialization of I_psi
			I_psi_1_1 = 0.0;
			I_psi_1_2 = 0.0;
			I_psi_1_3 = 0.0;
			I_psi_2_1 = 0.0;
			I_psi_2_2 = 0.0;
			I_psi_2_3 = 0.0;
			I_psi_3_1 = 0.0;
			I_psi_3_2 = 0.0;
			I_psi_3_3 = 0.0;
			I_psi_const = 0.0;
			//
			 for ( int n_psi = 0 ; n_psi <  N_psi ; n_psi++ )
			 {
				 // psi_A =< PSI <= psi_B
				 psi_1_A = 0.0;
				 psi_1_B = atan(Lq / Lp);
				 //
				 PSI_1 = ((psi_1_B-psi_1_A) / 2.0) * z_psi[n_psi] + (psi_1_B + psi_1_A) / 2.0;
				 //
				 J_psi_1 = (psi_1_B - psi_1_A) / 2.0;
                 //
				 psi_2_A = atan(Lq/Lp);
				 psi_2_B = M_PI / 2.0;
				 //
				 PSI_2 = ((psi_2_B - psi_2_A) / 2.0) * z_psi[n_psi] + (psi_2_B + psi_2_A) / 2.0;
				 //
				 J_psi_2 = (psi_2_B - psi_2_A) / 2.0;
				 //
				 B_1 = sqrt(a_v * pow(cos(PSI_1),2.0) + b_v * cos(PSI_1) * sin(PSI_1) + c_v * pow(sin(PSI_1),2.0));
				 //
				 B_2 = sqrt(a_v * pow(cos(PSI_2),2.0) + b_v * cos(PSI_2) * sin(PSI_2) + c_v * pow(sin(PSI_2),2.0));
				 //
				 Lp = (2.0 * sqrt(3.0)) / (sin(THETA_p) + sqrt(3.0) * cos(THETA_p));
				 L1 = Lp / cos(PSI_1);

				 Lq = (2.0 * sqrt(3.0)) / (sin(THETA_q) + sqrt(3.0) * cos(THETA_q));
				 L2 = Lq / sin(PSI_2);
				 //

				 K_functions(B_1, L1, K_1, ko);
				 K_functions(B_2, L2, K_2, ko);
                 //
				 Omega_f1_f1_1 = Omega_function_f1_f1(THETA_p, THETA_q, PSI_1, B_1, coef_f1_f1, K_1);
				 Omega_f1_f1_2 = Omega_function_f1_f1(THETA_p, THETA_q, PSI_2, B_2, coef_f1_f1, K_2);
				 //
				 Omega_f1_f2_1 = Omega_function_f1_f2(THETA_p, THETA_q, PSI_1, B_1, coef_f1_f2, K_1);
				 Omega_f1_f2_2 = Omega_function_f1_f2(THETA_p, THETA_q, PSI_2, B_2, coef_f1_f2, K_2);
				 //
				 Omega_f1_f3_1 = Omega_function_f1_f3(THETA_p, THETA_q, PSI_1, B_1, coef_f1_f3, K_1);
				 Omega_f1_f3_2 = Omega_function_f1_f3(THETA_p, THETA_q, PSI_2, B_2, coef_f1_f3, K_2);
				 //
				 Omega_f2_f1_1 = Omega_function_f2_f1(THETA_p, THETA_q, PSI_1, B_1, coef_f2_f1, K_1);
				 Omega_f2_f1_2 = Omega_function_f2_f1(THETA_p, THETA_q, PSI_2, B_2, coef_f2_f1, K_2);
				 //
				 Omega_f2_f2_1 = Omega_function_f2_f2(THETA_p, THETA_q, PSI_1, B_1, coef_f2_f2, K_1);
				 Omega_f2_f2_2 = Omega_function_f2_f2(THETA_p, THETA_q, PSI_2, B_2, coef_f2_f2, K_2);
				 //
				 Omega_f2_f3_1 = Omega_function_f2_f3(THETA_p, THETA_q, PSI_1, B_1, coef_f2_f3, K_1);
				 Omega_f2_f3_2 = Omega_function_f2_f3(THETA_p, THETA_q, PSI_2, B_2, coef_f2_f3, K_2);
				 //
				 Omega_f3_f1_1 = Omega_function_f3_f1(THETA_p, THETA_q, PSI_1, B_1, coef_f3_f1, K_1);
				 Omega_f3_f1_2 = Omega_function_f3_f1(THETA_p, THETA_q, PSI_2, B_2, coef_f3_f1, K_2);
				 //
				 Omega_f3_f2_1 = Omega_function_f3_f2(THETA_p, THETA_q, PSI_1, B_1, coef_f3_f2, K_1);
				 Omega_f3_f2_2 = Omega_function_f3_f2(THETA_p, THETA_q, PSI_2, B_2, coef_f3_f2, K_2);
				 //
				 Omega_f3_f3_1 = Omega_function_f3_f3(THETA_p, THETA_q, PSI_1, B_1, coef_f3_f3, K_1);
				 Omega_f3_f3_2 = Omega_function_f3_f3(THETA_p, THETA_q, PSI_2, B_2, coef_f3_f3, K_2);
				 //
				 Omega_const_1 = Omega_function_const(THETA_p, THETA_q, PSI_1, B_1, coef_const, K_1);
				 Omega_const_2 = Omega_function_const(THETA_p, THETA_q, PSI_2, B_2, coef_const, K_2);
				 //
         		 WPSI = w_psi[n_psi];
				 //
				 I_psi_1_1   +=  WPSI * (J_psi_1 * Omega_f1_f1_1 + J_psi_2 * Omega_f1_f1_2);
				 I_psi_1_2   +=  WPSI * (J_psi_1 * Omega_f1_f2_1 + J_psi_2 * Omega_f1_f2_2);
				 I_psi_1_3   +=  WPSI * (J_psi_1 * Omega_f1_f3_1 + J_psi_2 * Omega_f1_f3_2);
				 I_psi_2_1   +=  WPSI * (J_psi_1 * Omega_f2_f1_1 + J_psi_2 * Omega_f2_f1_2);
				 I_psi_2_2   +=  WPSI * (J_psi_1 * Omega_f2_f2_1 + J_psi_2 * Omega_f2_f2_2);
				 I_psi_2_3   +=  WPSI * (J_psi_1 * Omega_f2_f3_1 + J_psi_2 * Omega_f2_f3_2);
				 I_psi_3_1   +=  WPSI * (J_psi_1 * Omega_f3_f1_1 + J_psi_2 * Omega_f3_f1_2);
				 I_psi_3_2   +=  WPSI * (J_psi_1 * Omega_f3_f2_1 + J_psi_2 * Omega_f3_f2_2);
				 I_psi_3_3   +=  WPSI * (J_psi_1 * Omega_f3_f3_1 + J_psi_2 * Omega_f3_f3_2);
				 I_psi_const +=  WPSI * (J_psi_1 * Omega_const_1 + J_psi_2 * Omega_const_2);
			 }//end for ( int n_psi = 0 ; n_psi <  N_psi ; n_psi++ )
			 WTHETA_q = w_theta_q[n_theta_q];
			 //
			 I_theta_q_1_1     +=  WTHETA_q * I_psi_1_1;
			 I_theta_q_1_2     +=  WTHETA_q * I_psi_1_2;
			 I_theta_q_1_3     +=  WTHETA_q * I_psi_1_3;
			 I_theta_q_2_1     +=  WTHETA_q * I_psi_2_1;
			 I_theta_q_2_2     +=  WTHETA_q * I_psi_2_2;
			 I_theta_q_2_3     +=  WTHETA_q * I_psi_2_3;
			 I_theta_q_3_1     +=  WTHETA_q * I_psi_3_1;
			 I_theta_q_3_2     +=  WTHETA_q * I_psi_3_2;
			 I_theta_q_3_3     +=  WTHETA_q * I_psi_3_3;
			 I_theta_q_const   +=  WTHETA_q * I_psi_const;
			 
		 } //end for ( int n_theta_q = 0 ; n_theta_q <  N_theta_q ; n_theta_q++ )
		 //
		 I_theta_q_1_1     *=  J_theta_q ;
		 I_theta_q_1_2     *=  J_theta_q ;
		 I_theta_q_1_3     *=  J_theta_q ;
		 I_theta_q_2_1     *=  J_theta_q ;
		 I_theta_q_2_2     *=  J_theta_q ;
		 I_theta_q_2_3     *=  J_theta_q ;
		 I_theta_q_3_1     *=  J_theta_q ;
		 I_theta_q_3_2     *=  J_theta_q ;
		 I_theta_q_3_3     *=  J_theta_q ;
		 I_theta_q_const   *=  J_theta_q ;
		 //
		 WTHETA_p = w_theta_p[n_theta_p];
		 //
		 I_theta_p_1_1    += WTHETA_p * I_theta_q_1_1;
		 I_theta_p_1_2    += WTHETA_p * I_theta_q_1_2;
		 I_theta_p_1_3    += WTHETA_p * I_theta_q_1_3;
		 I_theta_p_2_1    += WTHETA_p * I_theta_q_2_1;
		 I_theta_p_2_2    += WTHETA_p * I_theta_q_2_2;
		 I_theta_p_2_3    += WTHETA_p * I_theta_q_2_3;
		 I_theta_p_3_1    += WTHETA_p * I_theta_q_3_1;
		 I_theta_p_3_2    += WTHETA_p * I_theta_q_3_2;
		 I_theta_p_3_3    += WTHETA_p * I_theta_q_3_3;
		 I_theta_p_const  += WTHETA_p * I_theta_q_const;

	 } //end for ( int n_theta_p = 0 ; n_theta_p <  N_theta_p ; n_theta_p++ )
	 //
	 I_theta_p_1_1     *=  J_theta_p ;
	 I_theta_p_1_2     *=  J_theta_p ;
	 I_theta_p_1_3     *=  J_theta_p ;
	 I_theta_p_2_1     *=  J_theta_p ;
	 I_theta_p_2_2     *=  J_theta_p ;
	 I_theta_p_2_3     *=  J_theta_p ;
	 I_theta_p_3_1     *=  J_theta_p ;
	 I_theta_p_3_2     *=  J_theta_p ;
	 I_theta_p_3_3     *=  J_theta_p ;
	 I_theta_p_const   *=  J_theta_p ;
	 //
	 double PRECOMP = 1.0 / 12.0;
	 // FINAL OUTPUT
	 I_DE[0] = PRECOMP * sqrt(vector_dot(r32,r32)) * sqrt(vector_dot(r54,r54)) * (Iunit * ko * I_theta_p_1_1 + 4.0 / (Iunit * ko) * I_theta_p_const);
	 I_DE[1] = PRECOMP * sqrt(vector_dot(r32,r32)) * sqrt(vector_dot(r51,r51)) * (Iunit * ko * I_theta_p_1_2 + 4.0 / (Iunit * ko) * I_theta_p_const);
	 I_DE[2] = PRECOMP * sqrt(vector_dot(r32,r32)) * sqrt(vector_dot(r41,r41)) * (Iunit * ko * I_theta_p_1_3 + 4.0 / (Iunit * ko) * I_theta_p_const);
	 I_DE[3] = PRECOMP * sqrt(vector_dot(r31,r31)) * sqrt(vector_dot(r54,r54)) * (Iunit * ko * I_theta_p_2_1 + 4.0 / (Iunit * ko) * I_theta_p_const);
	 I_DE[4] = PRECOMP * sqrt(vector_dot(r31,r31)) * sqrt(vector_dot(r51,r51)) * (Iunit * ko * I_theta_p_2_2 + 4.0 / (Iunit * ko) * I_theta_p_const);
	 I_DE[5] = PRECOMP * sqrt(vector_dot(r31,r31)) * sqrt(vector_dot(r41,r41)) * (Iunit * ko * I_theta_p_2_3 + 4.0 / (Iunit * ko) * I_theta_p_const);
	 I_DE[6] = PRECOMP * sqrt(vector_dot(r21,r21)) * sqrt(vector_dot(r54,r54)) * (Iunit * ko * I_theta_p_3_1 + 4.0 / (Iunit * ko) * I_theta_p_const);
	 I_DE[7] = PRECOMP * sqrt(vector_dot(r21,r21)) * sqrt(vector_dot(r51,r51)) * (Iunit * ko * I_theta_p_3_2 + 4.0 / (Iunit * ko) * I_theta_p_const);
	 I_DE[8] = PRECOMP * sqrt(vector_dot(r21,r21)) * sqrt(vector_dot(r41,r41)) * (Iunit * ko * I_theta_p_3_3 + 4.0 / (Iunit * ko) * I_theta_p_const);
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_g1_f1
// ***********************************************************************

void coefficients_const (const double r1[],const double r2[],const double r3[],const double r4[], const double r5[], const double ko, complex<double> coef[] )
{                         
	coef[0] = 1.0;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_g1_f1
// ***********************************************************************

void coefficients_f1_f1 (const double r1[],const double r2[],const double r3[],const double r4[], const double r5[], const double ko, complex<double> coef[] )
{                         
	double bb[3], cc[3], dd[3], ee[3];
	//
	for (int i = 0; i < 3; i++)
	{
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
			ee[i] = r5[i] - r1[i];
	}
	//
	complex<double> t1 = bb[0] * dd[0];
	complex<double> t3 = bb[1] * dd[1];
	complex<double> t5 = bb[2] * dd[2];
	complex<double> t7 = bb[0] * ee[0];
	complex<double> t9 = bb[1] * ee[1];
	complex<double> t11 = bb[2] * ee[2];
	complex<double> t13 = cc[0] * dd[0];
	complex<double> t15 = cc[1] * dd[1];
	complex<double> t17 = cc[2] * dd[2];
	complex<double> t25 = t1 / 0.12e2 + t3 / 0.12e2 + t5 / 0.12e2 - t7 / 0.6e1 - t9 / 0.6e1 - t11 / 0.6e1 - t13 / 0.6e1 - t15 / 0.6e1 - t17 / 0.6e1 + cc[0] * ee[0] / 0.3e1 + cc[1] * ee[1] / 0.3e1 + cc[2] * ee[2] / 0.3e1;
	complex<double> t27 = sqrt(0.3e1);
	complex<double> t29 = t1 * t27 / 0.12e2;
	complex<double> t31 = t3 * t27 / 0.12e2;
	complex<double> t35 = t5 * t27 / 0.12e2;
	//
	coef[0] = t25;
	coef[1] = t1 / 0.4e1 + t3 / 0.4e1 + t5 / 0.4e1;
	coef[2] = -t29 - t31 + t9 * t27 / 0.6e1 - t35 + t7 * t27 / 0.6e1 + t11 * t27 / 0.6e1;
	coef[3] = -t29 - t35 - t31 + t13 * t27 / 0.6e1 + t15 * t27 / 0.6e1 + t17 * t27 / 0.6e1;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f1_f2
// ***********************************************************************

void coefficients_f1_f2 (const double r1[],const double r2[],const double r3[],const double r4[], const double r5[], const double ko, complex<double> coef[] )
{                         
	double bb[3], cc[3], dd[3], ee[3];
	//
	for (int i = 0; i < 3; i++)
	{
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
			ee[i] = r5[i] - r1[i];
	}
	//
	complex<double> t1 = bb[1] * dd[1];
	complex<double> t3 = bb[0] * dd[0];
	complex<double> t5 = bb[2] * dd[2];
	complex<double> t7 = bb[1] * ee[1];
	complex<double> t9 = bb[2] * ee[2];
	complex<double> t11 = bb[0] * ee[0];
	complex<double> t25 = t1 / 0.12e2 + t3 / 0.12e2 + t5 / 0.12e2 - t7 / 0.6e1 - t9 / 0.6e1 - t11 / 0.6e1 - cc[0] * dd[0] / 0.6e1 - cc[1] * dd[1] / 0.6e1 - cc[2] * dd[2] / 0.6e1 + cc[1] * ee[1] / 0.3e1 + cc[2] * ee[2] / 0.3e1 + cc[0] * ee[0] / 0.3e1;
	complex<double> t26 = -t3 - t1 - t5;
	complex<double> t27 = sqrt(0.3e1);
	complex<double> t29 = t27 * bb[2] * dd[2];
	complex<double> t32 = t27 * bb[1] * dd[1];
	complex<double> t35 = t27 * cc[0] * dd[0];
	complex<double> t38 = t27 * bb[0] * dd[0];
	complex<double> t41 = t27 * cc[1] * dd[1];
	complex<double> t44 = t27 * cc[2] * dd[2];
	complex<double> t47 = t32 / 0.12e2;
	complex<double> t49 = t29 / 0.12e2;
	complex<double> t50 = t38 / 0.12e2;
	//
	coef[0] = t25;
	coef[1] = t26 / 0.2e1;
	coef[2] = t29 / 0.6e1 + t32 / 0.6e1 - t35 / 0.3e1 + t38 / 0.6e1 - t41 / 0.3e1 - t44 / 0.3e1;
	coef[3] = -t26 / 0.4e1;
	coef[4] = -t47 + t41 / 0.6e1 - t49 - t50 + t44 / 0.6e1 + t35 / 0.6e1;
	coef[5] = -t50 - t47 + t9 * t27 / 0.6e1 - t49 + t11 * t27 / 0.6e1 + t7 * t27 / 0.6e1;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f1_f3
// ***********************************************************************

void coefficients_f1_f3 (const double r1[],const double r2[],const double r3[],const double r4[], const double r5[], const double ko, complex<double> coef[] )
{                         
	double bb[3], cc[3], dd[3], ee[3];
	//
	for (int i = 0; i < 3; i++)
	{
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
			ee[i] = r5[i] - r1[i];
	}
	//
	complex<double> t1 = bb[1] * dd[1];
	complex<double> t3 = bb[0] * dd[0];
	complex<double> t5 = bb[2] * dd[2];
	complex<double> t7 = bb[1] * ee[1];
	complex<double> t9 = bb[2] * ee[2];
	complex<double> t11 = bb[0] * ee[0];
	complex<double> t13 = cc[0] * dd[0];
	complex<double> t15 = cc[1] * dd[1];
	complex<double> t17 = cc[2] * dd[2];
	complex<double> t25 = t1 / 0.12e2 + t3 / 0.12e2 + t5 / 0.12e2 - t7 / 0.6e1 - t9 / 0.6e1 - t11 / 0.6e1 - t13 / 0.6e1 - t15 / 0.6e1 - t17 / 0.6e1 + cc[1] * ee[1] / 0.3e1 + cc[2] * ee[2] / 0.3e1 + cc[0] * ee[0] / 0.3e1;
	complex<double> t27 = sqrt(0.3e1);
	complex<double> t30 = t27 * bb[2] * ee[2] / 0.6e1;
	complex<double> t33 = t27 * bb[1] * ee[1] / 0.6e1;
	complex<double> t39 = t27 * bb[0] * ee[0] / 0.6e1;
	complex<double> t49 = t1 * t27 / 0.12e2;
	complex<double> t51 = t5 * t27 / 0.12e2;
	complex<double> t55 = t3 * t27 / 0.12e2;
	//
	coef[0] = t25;
	coef[1] = -t11 / 0.2e1 - t7 / 0.2e1 - t9 / 0.2e1;
	coef[2] = t30 + t33 - t27 * cc[0] * ee[0] / 0.3e1 + t39 - t27 * cc[1] * ee[1] / 0.3e1 - t27 * cc[2] * ee[2] / 0.3e1;
	coef[3] = t1 / 0.4e1 + t3 / 0.4e1 + t5 / 0.4e1;
	coef[4] = -t49 - t51 + t15 * t27 / 0.6e1 - t55 + t13 * t27 / 0.6e1 + t17 * t27 / 0.6e1;
	coef[5] = -t51 - t49 + t39 - t55 + t33 + t30;

}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f2_f2
// ***********************************************************************

void coefficients_f2_f1 (const double r1[],const double r2[],const double r3[],const double r4[], const double r5[], const double ko, complex<double> coef[] )
{                         
	double aa[3], bb[3], cc[3], dd[3], ee[3];
	//
	for (int i = 0; i < 3; i++)
	{
			aa[i] = r1[i] - r1[i];
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
			ee[i] = r5[i] - r1[i];
	}
	//complex<double> j   = Iunit;
	//
	complex<double> t3 = bb[0] * dd[0];
	complex<double> t7 = bb[1] * dd[1];
	complex<double> t11 = bb[2] * dd[2];
	complex<double> t13 = cc[1] * dd[1];
	complex<double> t19 = cc[0] * dd[0];
	complex<double> t23 = cc[2] * dd[2];
	complex<double> t25 = cc[0] * ee[0] / 0.3e1 + t3 / 0.12e2 - bb[0] * ee[0] / 0.6e1 + t7 / 0.12e2 - bb[2] * ee[2] / 0.6e1 + t11 / 0.12e2 - t13 / 0.6e1 - bb[1] * ee[1] / 0.6e1 + cc[2] * ee[2] / 0.3e1 - t19 / 0.6e1 + cc[1] * ee[1] / 0.3e1 - t23 / 0.6e1;
	complex<double> t26 = -t3 - t7 - t11;
	complex<double> t27 = sqrt(0.3e1);
	complex<double> t28 = t27 * bb[1];
	complex<double> t29 = t28 * dd[1];
	complex<double> t30 = t29 / 0.12e2;
	complex<double> t31 = t27 * bb[0];
	complex<double> t32 = t31 * ee[0];
	complex<double> t34 = t31 * dd[0];
	complex<double> t35 = t34 / 0.12e2;
	complex<double> t36 = t27 * bb[2];
	complex<double> t37 = t36 * ee[2];
	complex<double> t39 = t28 * ee[1];
	complex<double> t41 = t36 * dd[2];
	complex<double> t42 = t41 / 0.12e2;
	//
	coef[0] = t25;
	coef[1] = t26 / 0.2e1;
	coef[2] = -t26 / 0.4e1;
	coef[3] = -t30 + t32 / 0.6e1 - t35 + t37 / 0.6e1 + t39 / 0.6e1 - t42;
	coef[4] = t23 * t27 / 0.6e1 - t30 + t19 * t27 / 0.6e1 - t42 + t13 * t27 / 0.6e1 - t35;
	coef[5] = t29 / 0.6e1 + t34 / 0.6e1 - t39 / 0.3e1 + t41 / 0.6e1 - t32 / 0.3e1 - t37 / 0.3e1;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f2_f2
// ***********************************************************************

void coefficients_f2_f2 (const double r1[],const double r2[],const double r3[],const double r4[], const double r5[], const double ko, complex<double> coef[] )
{                         
	double bb[3], cc[3], dd[3], ee[3];
	//
	for (int i = 0; i < 3; i++)
	{
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
			ee[i] = r5[i] - r1[i];
	}
	//
	complex<double> t1 = bb[2] * dd[2];
	complex<double> t2 = bb[0] * dd[0];
	complex<double> t3 = bb[1] * dd[1];
	complex<double> t4 = t1 + t2 + t3;
	complex<double> t26 = -bb[1] * ee[1] / 0.6e1 + t3 / 0.12e2 + t1 / 0.12e2 - bb[2] * ee[2] / 0.6e1 + t2 / 0.12e2 - bb[0] * ee[0] / 0.6e1 - cc[0] * dd[0] / 0.6e1 - cc[1] * dd[1] / 0.6e1 - cc[2] * dd[2] / 0.6e1 + cc[0] * ee[0] / 0.3e1 + cc[1] * ee[1] / 0.3e1 + cc[2] * ee[2] / 0.3e1;
	complex<double> t27 = sqrt(0.3e1);
	complex<double> t28 = t27 * bb[0];
	complex<double> t29 = t28 * dd[0];
	complex<double> t30 = t29 / 0.6e1;
	complex<double> t31 = t27 * bb[1];
	complex<double> t32 = t31 * dd[1];
	complex<double> t33 = t32 / 0.6e1;
	complex<double> t35 = t27 * cc[0] * dd[0];
	complex<double> t38 = t27 * cc[1] * dd[1];
	complex<double> t40 = t27 * bb[2];
	complex<double> t41 = t40 * dd[2];
	complex<double> t42 = t41 / 0.6e1;
	complex<double> t44 = t27 * cc[2] * dd[2];
	complex<double> t47 = t40 * ee[2];
	complex<double> t49 = t28 * ee[0];
	complex<double> t51 = t31 * ee[1];
	complex<double> t55 = t41 / 0.12e2;
	complex<double> t57 = t29 / 0.12e2;
	complex<double> t58 = t32 / 0.12e2;
	//
	coef[0] = t4;
	coef[1] = t26;
	coef[2] = -t4 / 0.2e1;
	coef[3] = t30 + t33 - t35 / 0.3e1 - t38 / 0.3e1 + t42 - t44 / 0.3e1;
	coef[4] = -t4 / 0.2e1;
	coef[5] = -t47 / 0.3e1 + t30 + t33 + t42 - t49 / 0.3e1 - t51 / 0.3e1;
	coef[6] = t51 / 0.6e1 - t55 + t47 / 0.6e1 - t57 - t58 + t49 / 0.6e1;
	coef[7] = t35 / 0.6e1 + t44 / 0.6e1 - t55 - t57 - t58 + t38 / 0.6e1;
	coef[8] = t4 / 0.4e1;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f2_f3
// ***********************************************************************

void coefficients_f2_f3 (const double r1[],const double r2[],const double r3[],const double r4[], const double r5[], const double ko, complex<double> coef[] )
{                         
	double bb[3], cc[3], dd[3], ee[3];
	//
	for (int i = 0; i < 3; i++)
	{
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
			ee[i] = r5[i] - r1[i];
	}
	//
	complex<double> t1 = bb[2] * ee[2];
	complex<double> t2 = bb[0] * ee[0];
	complex<double> t3 = bb[1] * ee[1];
	complex<double> t4 = t1 + t2 + t3;
	complex<double> t6 = bb[2] * dd[2];
	complex<double> t8 = bb[1] * dd[1];
	complex<double> t11 = bb[0] * dd[0];
	complex<double> t14 = cc[0] * dd[0];
	complex<double> t16 = cc[1] * dd[1];
	complex<double> t18 = cc[2] * dd[2];
	complex<double> t26 = -t1 / 0.6e1 + t6 / 0.12e2 + t8 / 0.12e2 - t2 / 0.6e1 + t11 / 0.12e2 - t3 / 0.6e1 - t14 / 0.6e1 - t16 / 0.6e1 - t18 / 0.6e1 + cc[0] * ee[0] / 0.3e1 + cc[1] * ee[1] / 0.3e1 + cc[2] * ee[2] / 0.3e1;
	complex<double> t27 = t8 + t11 + t6;
	complex<double> t28 = sqrt(0.3e1);
	complex<double> t31 = t28 * bb[2];
	complex<double> t32 = t31 * dd[2];
	complex<double> t33 = t32 / 0.12e2;
	complex<double> t34 = t28 * bb[1];
	complex<double> t35 = t34 * dd[1];
	complex<double> t36 = t35 / 0.12e2;
	complex<double> t41 = t28 * bb[0];
	complex<double> t42 = t41 * dd[0];
	complex<double> t43 = t42 / 0.12e2;
	complex<double> t45 = t31 * ee[2];
	complex<double> t46 = t45 / 0.6e1;
	complex<double> t47 = t34 * ee[1];
	complex<double> t48 = t47 / 0.6e1;
	complex<double> t49 = t41 * ee[0];
	complex<double> t50 = t49 / 0.6e1;
	//
	coef[0] = t4;
	coef[1] = t26;
	coef[2] = t27 / 0.4e1;
	coef[3] = t18 * t28 / 0.6e1 - t33 - t36 + t14 * t28 / 0.6e1 + t16 * t28 / 0.6e1 - t43;
	coef[4] = -t43 + t46 + t48 - t36 - t33 + t50;
	coef[5] = t46 + t50 - t28 * cc[0] * ee[0] / 0.3e1 - t28 * cc[1] * ee[1] / 0.3e1 + t48 - t28 * cc[2] * ee[2] / 0.3e1;
	coef[6] = -t27 / 0.2e1;
	coef[7] = -t47 / 0.3e1 + t32 / 0.6e1 + t35 / 0.6e1 + t42 / 0.6e1 - t45 / 0.3e1 - t49 / 0.3e1;
	coef[8] = -t4 / 0.2e1;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f3_f1
// ***********************************************************************

void coefficients_f3_f1 (const double r1[],const double r2[],const double r3[],const double r4[], const double r5[], const double ko, complex<double> coef[] )
{                         
	double bb[3], cc[3], dd[3], ee[3];
	//
	for (int i = 0; i < 3; i++)
	{
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
			ee[i] = r5[i] - r1[i];
	}
	//
	complex<double> t1 = bb[0] * dd[0];
	complex<double> t3 = bb[1] * dd[1];
	complex<double> t5 = bb[2] * dd[2];
	complex<double> t7 = bb[0] * ee[0];
	complex<double> t9 = bb[1] * ee[1];
	complex<double> t11 = bb[2] * ee[2];
	complex<double> t13 = cc[0] * dd[0];
	complex<double> t15 = cc[1] * dd[1];
	complex<double> t17 = cc[2] * dd[2];
	complex<double> t25 = t1 / 0.12e2 + t3 / 0.12e2 + t5 / 0.12e2 - t7 / 0.6e1 - t9 / 0.6e1 - t11 / 0.6e1 - t13 / 0.6e1 - t15 / 0.6e1 - t17 / 0.6e1 + cc[2] * ee[2] / 0.3e1 + cc[0] * ee[0] / 0.3e1 + cc[1] * ee[1] / 0.3e1;
	complex<double> t26 = sqrt(0.3e1);
	complex<double> t27 = t26 * cc[1];
	complex<double> t29 = t27 * dd[1] / 0.6e1;
	complex<double> t30 = t26 * cc[0];
	complex<double> t32 = t30 * dd[0] / 0.6e1;
	complex<double> t35 = t26 * cc[2];
	complex<double> t37 = t35 * dd[2] / 0.6e1;
	complex<double> t46 = t1 * t26 / 0.12e2;
	complex<double> t48 = t3 * t26 / 0.12e2;
	complex<double> t52 = t5 * t26 / 0.12e2;
	//
	coef[0] = t25;
	coef[1] = t29 + t32 - t27 * ee[1] / 0.3e1 + t37 - t30 * ee[0] / 0.3e1 - t35 * ee[2] / 0.3e1;
	coef[2] = -t17 / 0.2e1 - t13 / 0.2e1 - t15 / 0.2e1;
	coef[3] = t1 / 0.4e1 + t3 / 0.4e1 + t5 / 0.4e1;
	coef[4] = -t46 - t48 + t9 * t26 / 0.6e1 - t52 + t7 * t26 / 0.6e1 + t11 * t26 / 0.6e1;
	coef[5] = -t46 - t52 - t48 + t32 + t29 + t37;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f3_f2
// ***********************************************************************

void coefficients_f3_f2 (const double r1[],const double r2[],const double r3[],const double r4[], const double r5[], const double ko, complex<double> coef[] )
{                         
	double bb[3], cc[3], dd[3], ee[3];
	//
	for (int i = 0; i < 3; i++)
	{
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
			ee[i] = r5[i] - r1[i];
	}
	//
	complex<double> t1 = cc[2] * dd[2];
	complex<double> t2 = cc[0] * dd[0];
	complex<double> t3 = cc[1] * dd[1];
	complex<double> t4 = t1 + t2 + t3;
	complex<double> t5 = bb[2] * dd[2];
	complex<double> t7 = bb[1] * ee[1];
	complex<double> t9 = bb[2] * ee[2];
	complex<double> t11 = bb[0] * dd[0];
	complex<double> t15 = bb[1] * dd[1];
	complex<double> t19 = bb[0] * ee[0];
	complex<double> t26 = t5 / 0.12e2 - t7 / 0.6e1 - t9 / 0.6e1 + t11 / 0.12e2 + cc[1] * ee[1] / 0.3e1 + t15 / 0.12e2 - t2 / 0.6e1 - t3 / 0.6e1 - t19 / 0.6e1 - t1 / 0.6e1 + cc[2] * ee[2] / 0.3e1 + cc[0] * ee[0] / 0.3e1;
	complex<double> t27 = t15 + t11 + t5;
	complex<double> t28 = sqrt(0.3e1);
	complex<double> t30 = t28 * bb[1] * dd[1];
	complex<double> t31 = t30 / 0.12e2;
	complex<double> t35 = t28 * bb[2] * dd[2];
	complex<double> t36 = t35 / 0.12e2;
	complex<double> t42 = t28 * bb[0] * dd[0];
	complex<double> t43 = t42 / 0.12e2;
	complex<double> t45 = t28 * cc[1];
	complex<double> t46 = t45 * dd[1];
	complex<double> t47 = t46 / 0.6e1;
	complex<double> t48 = t28 * cc[0];
	complex<double> t49 = t48 * dd[0];
	complex<double> t50 = t49 / 0.6e1;
	complex<double> t51 = t28 * cc[2];
	complex<double> t52 = t51 * dd[2];
	complex<double> t53 = t52 / 0.6e1;
	//
	coef[0] = t4;
	coef[1] = t26;
	coef[2] = -t4 / 0.2e1;
	coef[3] = t27 / 0.4e1;
	coef[4] = -t31 + t9 * t28 / 0.6e1 - t36 + t7 * t28 / 0.6e1 + t19 * t28 / 0.6e1 - t43;
	coef[5] = -t43 + t47 + t50 - t36 + t53 - t31;
	coef[6] = -t48 * ee[0] / 0.3e1 + t47 + t50 + t53 - t51 * ee[2] / 0.3e1 - t45 * ee[1] / 0.3e1;
	coef[7] = -t27 / 0.2e1;
	coef[8] = t35 / 0.6e1 - t49 / 0.3e1 + t30 / 0.6e1 + t42 / 0.6e1 - t52 / 0.3e1 - t46 / 0.3e1;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f3_f3
// ***********************************************************************

void coefficients_f3_f3 (const double r1[],const double r2[],const double r3[],const double r4[], const double r5[], const double ko, complex<double> coef[] )
{                         
	double bb[3], cc[3], dd[3], ee[3];
	//
	for (int i = 0; i < 3; i++)
	{
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
			ee[i] = r5[i] - r1[i];
	}
	//
	complex<double> t1 = cc[2] * ee[2];
	complex<double> t2 = cc[0] * ee[0];
	complex<double> t3 = cc[1] * ee[1];
	complex<double> t5 = bb[1] * ee[1];
	complex<double> t7 = bb[2] * dd[2];
	complex<double> t9 = bb[1] * dd[1];
	complex<double> t11 = bb[2] * ee[2];
	complex<double> t13 = bb[0] * dd[0];
	complex<double> t15 = bb[0] * ee[0];
	complex<double> t17 = cc[0] * dd[0];
	complex<double> t19 = cc[1] * dd[1];
	complex<double> t21 = cc[2] * dd[2];
	complex<double> t26 = -t5 / 0.6e1 + t7 / 0.12e2 + t9 / 0.12e2 - t11 / 0.6e1 + t13 / 0.12e2 - t15 / 0.6e1 - t17 / 0.6e1 - t19 / 0.6e1 - t21 / 0.6e1 + t1 / 0.3e1 + t2 / 0.3e1 + t3 / 0.3e1;
	complex<double> t27 = sqrt(0.3e1);
	complex<double> t28 = t27 * cc[2];
	complex<double> t30 = t28 * ee[2] / 0.3e1;
	complex<double> t32 = t28 * dd[2] / 0.6e1;
	complex<double> t33 = t27 * cc[0];
	complex<double> t35 = t33 * dd[0] / 0.6e1;
	complex<double> t36 = t27 * cc[1];
	complex<double> t38 = t36 * dd[1] / 0.6e1;
	complex<double> t40 = t33 * ee[0] / 0.3e1;
	complex<double> t42 = t36 * ee[1] / 0.3e1;
	complex<double> t46 = t27 * bb[1] * ee[1] / 0.6e1;
	complex<double> t49 = t27 * bb[2] * ee[2] / 0.6e1;
	complex<double> t52 = t27 * bb[0] * ee[0] / 0.6e1;
	complex<double> t56 = t9 * t27 / 0.12e2;
	complex<double> t58 = t7 * t27 / 0.12e2;
	complex<double> t60 = t13 * t27 / 0.12e2;
	//
	coef[0] = t1 + t2 + t3;
	coef[1] = t26;
	coef[2] = -t30 + t32 + t35 + t38 - t40 - t42;
	coef[3] = t46 + t49 - t42 - t30 + t52 - t40;
	coef[4] = t9 / 0.4e1 + t13 / 0.4e1 + t7 / 0.4e1;
	coef[5] = -t56 + t49 - t58 - t60 + t46 + t52;
	coef[6] = -t56 - t60 + t38 - t58 + t35 + t32;
	coef[7] = -t5 / 0.2e1 - t11 / 0.2e1 - t15 / 0.2e1;
	coef[8] = -t21 / 0.2e1 - t17 / 0.2e1 - t19 / 0.2e1;
}

// ***********************************************************************
//			IMPLEMENTATION OF void K_functions
// ***********************************************************************

void K_functions (double GAMMA, double L, complex<double> K[], const double ko)
{
	complex<double> a  = Iunit * ko * GAMMA;
	//
	K[0] = (1.0 / pow(a,2.0) ) * (1.0 - exp(-a * L) - a * L * exp(-a * L) );
	K[1] = (1.0 / pow(a,3.0) ) * (2.0 - 2.0 * exp(-a * L) - 2.0 * a * L * exp(-a * L) - pow(a,2.0) * pow(L,2.0) * exp(-a * L) );
	K[2] = -(pow(L,3.0) * exp(-a * L) ) / a + (3.0 / a) * K[1];
	K[3] = -(pow(L,4.0) * exp(-a * L) ) / a + (4.0 / a) * K[2];
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> Omega_function_const
// ***********************************************************************

complex<double> Omega_function_const (double theta_p, double theta_q, double Psi, double GAMMA, complex<double> coef[], complex<double> K[])
{                                 
	complex<double> X;
	//
	X = K[1] * (sin(Psi) * cos(Psi) * coef[0]) / GAMMA;
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> Omega_function_f1_f1
// ***********************************************************************


complex<double> Omega_function_f1_f1 (double theta_p, double theta_q, double Psi, double GAMMA, complex<double> coef[], complex<double> K[])
{                                 
	complex<double> X;
	//
	double t2 = sin(Psi);
	double t3 = t2 * t2;
	double t4 = cos(Psi);
	double t5 = t4 * t4;
	double t8 = t3 * t5 / GAMMA;
	double t9 = sin(theta_p);
	double t10 = sin(theta_q);
	double t15 = cos(theta_q);
	double t16 = cos(theta_p);
	//
	X = K[3] * (t8 * t9 * t10 * coef[0] + t8 * t15 * t16 * coef[1] + t8 * t10 * t16 * coef[2] + t8 * t15 * t9 * coef[3]);
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> Omega_function_f1_f2
// ***********************************************************************

complex<double> Omega_function_f1_f2 (double theta_p, double theta_q, double Psi, double GAMMA, complex<double> coef[], complex<double> K[])
{                                 
	complex<double> X;
	//
	complex<double> t2 = sin(Psi);
	complex<double> t3 = cos(Psi);
	complex<double> t4 = t3 * t3;
	complex<double> t5 = t2 * t4;
	complex<double> t6 = 0.1e1 / GAMMA;
	complex<double> t7 = cos(theta_p);
	complex<double> t12 = sin(theta_p);
	complex<double> t20 = t2 * t2;
	complex<double> t22 = t20 * t4 * t6;
	complex<double> t23 = sin(theta_q);
	complex<double> t28 = cos(theta_q);
	//
	X = K[2] * (t5 * t6 * t7 * coef[1] + t5 * t6 * t12 * coef[2]) + K[3] * (t22 * t12 * t23 * coef[0] + t22 * t28 * t7 * coef[3] + t22 * t28 * t12 * coef[4] + t22 * t23 * t7 * coef[5]);
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> Omega_function_f1_f3
// ***********************************************************************

complex<double> Omega_function_f1_f3 (double theta_p, double theta_q, double Psi, double GAMMA, complex<double> coef[], complex<double> K[])
{                                 
	complex<double> X;
	//
	complex<double> t2 = sin(Psi);
	complex<double> t3 = cos(Psi);
	complex<double> t4 = t3 * t3;
	complex<double> t5 = t2 * t4;
	complex<double> t6 = 0.1e1 / GAMMA;
	complex<double> t7 = cos(theta_p);
	complex<double> t12 = sin(theta_p);
	complex<double> t20 = t2 * t2;
	complex<double> t22 = t20 * t4 * t6;
	complex<double> t23 = sin(theta_q);
	complex<double> t28 = cos(theta_q);
	//
	X = K[2] * (t5 * t6 * t7 * coef[1] + t5 * t6 * t12 * coef[2]) + K[3] * (t22 * t12 * t23 * coef[0] + t22 * t28 * t7 * coef[3] + t22 * t28 * t12 * coef[4] + t22 * t23 * t7 * coef[5]);
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> Omega_function_f2_f1
// ***********************************************************************

complex<double> Omega_function_f2_f1 (double theta_p, double theta_q, double Psi, double GAMMA, complex<double> coef[], complex<double> K[])
{                                 
	complex<double> X;
	//
	complex<double> t2 = sin(Psi);
	complex<double> t3 = t2 * t2;
	complex<double> t4 = cos(Psi);
	complex<double> t5 = t3 * t4;
	complex<double> t6 = 0.1e1 / GAMMA;
	complex<double> t7 = cos(theta_q);
	complex<double> t12 = sin(theta_q);
	complex<double> t20 = t4 * t4;
	complex<double> t22 = t3 * t20 * t6;
	complex<double> t23 = sin(theta_p);
	complex<double> t28 = cos(theta_p);
	//
	X = K[2] * (t5 * t6 * t7 * coef[1] + t5 * t6 * t12 * coef[5]) + K[3] * (t22 * t23 * t12 * coef[0] + t22 * t7 * t28 * coef[2] + t22 * t12 * t28 * coef[3] + t22 * t7 * t23 * coef[4]);
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> Omega_function_f2_f2
// ***********************************************************************

complex<double> Omega_function_f2_f2 (double theta_p, double theta_q, double Psi, double GAMMA, complex<double> coef[], complex<double> K[])
{                                 
	complex<double> X;
	//
	complex<double> t2 = sin(Psi);
	complex<double> t3 = t2 * t2;
	complex<double> t4 = cos(Psi);
	complex<double> t5 = t3 * t4;
	complex<double> t6 = 0.1e1 / GAMMA;
	complex<double> t7 = cos(theta_q);
	complex<double> t12 = sin(theta_q);
	complex<double> t17 = t4 * t4;
	complex<double> t18 = t2 * t17;
	complex<double> t19 = cos(theta_p);
	complex<double> t24 = sin(theta_p);
	complex<double> t39 = t3 * t17 * t6;
	//
	X = K[2] * (t5 * t6 * t7 * coef[4] + t5 * t6 * t12 * coef[5] + t18 * t6 * t19 * coef[2] + t18 * t6 * t24 * coef[3]) + K[1] * t2 * t4 * t6 * coef[0] + K[3] * (t39 * t24 * t12 * coef[1] + t39 * t19 * t12 * coef[6] + t39 * t7 * t24 * coef[7] + t39 * t7 * t19 * coef[8]);
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> Omega_function_f2_f3
// ***********************************************************************

complex<double> Omega_function_f2_f3 (double theta_p, double theta_q, double Psi, double GAMMA, complex<double> coef[], complex<double> K[])
{                                 
	complex<double> X;
	//
	complex<double> t2 = sin(Psi);
	complex<double> t3 = cos(Psi);
	complex<double> t4 = t3 * t3;
	complex<double> t5 = t2 * t4;
	complex<double> t6 = 0.1e1 / GAMMA;
	complex<double> t7 = cos(theta_p);
	complex<double> t12 = sin(theta_p);
	complex<double> t17 = t2 * t2;
	complex<double> t18 = t17 * t3;
	complex<double> t19 = cos(theta_q);
	complex<double> t24 = sin(theta_q);
	complex<double> t39 = t17 * t4 * t6;
	//
	X = K[2] * (t5 * t6 * t7 * coef[8] + t5 * t6 * t12 * coef[5] + t18 * t6 * t19 * coef[6] + t18 * t6 * t24 * coef[7]) + K[1] * t2 * t3 * t6 * coef[0] + K[3] * (t39 * t12 * t24 * coef[1] + t39 * t19 * t7 * coef[2] + t39 * t19 * t12 * coef[3] + t39 * t24 * t7 * coef[4]);
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> Omega_function_f3_f1
// ***********************************************************************

complex<double> Omega_function_f3_f1 (double theta_p, double theta_q, double Psi, double GAMMA, complex<double> coef[], complex<double> K[])
{                                 
	complex<double> X;
	//
	complex<double> t2 = sin(Psi);
	complex<double> t3 = t2 * t2;
	complex<double> t4 = cos(Psi);
	complex<double> t5 = t3 * t4;
	complex<double> t6 = 0.1e1 / GAMMA;
	complex<double> t7 = sin(theta_q);
	complex<double> t12 = cos(theta_q);
	complex<double> t20 = t4 * t4;
	complex<double> t22 = t3 * t20 * t6;
	complex<double> t23 = sin(theta_p);
	complex<double> t28 = cos(theta_p);
	//
	X = K[2] * (t5 * t6 * t7 * coef[1] + t5 * t6 * t12 * coef[2]) + K[3] * (t22 * t23 * t7 * coef[0] + t22 * t12 * t28 * coef[3] + t22 * t7 * t28 * coef[4] + t22 * t12 * t23 * coef[5]);
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> Omega_function_f3_f2
// ***********************************************************************

complex<double> Omega_function_f3_f2 (double theta_p, double theta_q, double Psi, double GAMMA, complex<double> coef[], complex<double> K[])
{                                 
	complex<double> X;
	//
	complex<double> t2 = sin(Psi);
	complex<double> t3 = cos(Psi);
	complex<double> t4 = t3 * t3;
	complex<double> t5 = t2 * t4;
	complex<double> t6 = 0.1e1 / GAMMA;
	complex<double> t7 = sin(theta_p);
	complex<double> t12 = t2 * t2;
	complex<double> t13 = t12 * t3;
	complex<double> t14 = cos(theta_q);
	complex<double> t19 = cos(theta_p);
	complex<double> t24 = sin(theta_q);
	complex<double> t39 = t12 * t4 * t6;
	//
	X = K[2] * (t5 * t6 * t7 * coef[8] + t13 * t6 * t14 * coef[2] + t5 * t6 * t19 * coef[7] + t13 * t6 * t24 * coef[6]) + K[1] * t2 * t3 * t6 * coef[0] + K[3] * (t39 * t7 * t24 * coef[1] + t39 * t14 * t7 * coef[5] + t39 * t14 * t19 * coef[3] + t39 * t24 * t19 * coef[4]);
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> Omega_function_f3_f3
// ***********************************************************************

complex<double> Omega_function_f3_f3 (double theta_p, double theta_q, double Psi, double GAMMA, complex<double> coef[], complex<double> K[])
{                                 
	complex<double> X;
	//
	complex<double> t2 = sin(Psi);
	complex<double> t3 = t2 * t2;
	complex<double> t4 = cos(Psi);
	complex<double> t5 = t3 * t4;
	complex<double> t6 = 0.1e1 / GAMMA;
	complex<double> t7 = cos(theta_q);
	complex<double> t12 = sin(theta_q);
	complex<double> t17 = t4 * t4;
	complex<double> t18 = t2 * t17;
	complex<double> t19 = sin(theta_p);
	complex<double> t24 = cos(theta_p);
	complex<double> t39 = t3 * t17 * t6;
	//
	X = K[2] * (t5 * t6 * t7 * coef[8] + t5 * t6 * t12 * coef[2] + t18 * t6 * t19 * coef[3] + t18 * t6 * t24 * coef[7]) + K[1] * t2 * t4 * t6 * coef[0] + K[3] * (t39 * t19 * t12 * coef[1] + t39 * t12 * t24 * coef[5] + t39 * t7 * t19 * coef[6] + t39 * t7 * t24 * coef[4]);
	// Final Output
	return X;
}
