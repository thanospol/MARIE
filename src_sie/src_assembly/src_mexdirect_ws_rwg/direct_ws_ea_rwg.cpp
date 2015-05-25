/**************************************************************************************************************************
           
		           DIRECT_WS_EA_RWG.cpp

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
  r1,r2,r3,r4 = point vectors of the triangular element's vertices
  Outer triangle P:(rp1,rp2,rp3)=(r1,r2,r3)
  Inner triangle Q:(rq1,rq2,rq3)=(r2,r1,r4)
  N_theta,N_psi = order of the Gauss-Legendre quadrature for both dimensions of the remaining 2-D smooth integral

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

#include "direct_ws_ea_rwg.h"

// ***********************************************************************
//			IMPLEMENTATION OF void DIRECT
// ***********************************************************************

void direct_ws_ea_rwg (const double r1[], const double r2[] ,const double r3[], const double r4[], const double ko, const int N_theta, const int N_psi, const double w_theta[], const double z_theta[], const double w_psi[], const double z_psi[], complex<double> I_DE[] )
{
	// ************************************************
	//			DECLARATION OF KEY VARIABLES
	// ************************************************

	// 1. Various

	double alpha[3], beta[3], gamma[3], r21[3], r31[3], r32[3], r41[3], r42[3];
	//
	complex<double> Ising_const, Ising_f1_f1, Ising_f1_f2, Ising_f1_f3, Ising_f2_f1, Ising_f2_f2, Ising_f2_f3, Ising_f3_f1, Ising_f3_f2, Ising_f3_f3;
	//
	double THETA;
	double J_theta;
	double tPsiA;
	double tPsiB;
	double PsiA;
	double PsiB;
	double psi_A;
	double psi_B;
	//
	complex<double> N[12], Nm[12];
	//
	double PSI;
	double J_psi;
	double B;
	double Bm;
	//
	double b0;
	double b1;
	double b2;
	double B1;
	double B2;
	//
	double WPSI, WTHETA;
	//
	complex<double> I_theta_const, I_theta_f1_f1, I_theta_f1_f2, I_theta_f1_f3, I_theta_f2_f1, I_theta_f2_f2, I_theta_f2_f3, I_theta_f3_f1, I_theta_f3_f2, I_theta_f3_f3;
	complex<double> I_psi_const, I_psi_f1_f1, I_psi_f1_f2, I_psi_f1_f3, I_psi_f2_f1, I_psi_f2_f2, I_psi_f2_f3, I_psi_f3_f1, I_psi_f3_f2, I_psi_f3_f3;
	//
	complex<double> X_const, X_f1_f1, X_f1_f2, X_f1_f3, X_f2_f1, X_f2_f2, X_f2_f3, X_f3_f1, X_f3_f2, X_f3_f3;
	//
	complex<double> I_const[6], I_f1_f1[6], I_f1_f2[6], I_f1_f3[6], I_f2_f1[6], I_f2_f2[6], I_f2_f3[6], I_f3_f1[6], I_f3_f2[6], I_f3_f3[6]; 
	//
	complex<double> Iconst, If1_f1, If1_f2, If1_f3, If2_f1, If2_f2, If2_f3, If3_f1, If3_f2, If3_f3; 

	// 2. Coefficients' parameters

	complex<double> coef_const[3],  coefm_const[3];
	complex<double> coef_f1_f1[40], coefm_f1_f1[40];
	complex<double> coef_f1_f2[43], coefm_f1_f2[43];
	complex<double> coef_f1_f3[43], coefm_f1_f3[43];
	complex<double> coef_f2_f1[43], coefm_f2_f1[43];
	complex<double> coef_f2_f2[40], coefm_f2_f2[40];
	complex<double> coef_f2_f3[43], coefm_f2_f3[43];
	complex<double> coef_f3_f1[43], coefm_f3_f1[43];
	complex<double> coef_f3_f2[43], coefm_f3_f2[43];
	complex<double> coef_f3_f3[43], coefm_f3_f3[43];
	// 4.

	double theta_A;
	double theta_B;

	// ************************************************
	//			MAIN CODE
	// ************************************************
	// Evaluate alpha, beta, gamma r21, r31, r32, r41, r42 parameters
	for (int i = 0; i < 3; i++)
	{
			alpha[i] = (r2[i] - r1[i]) / 2;
			beta[i]  = (2 * r3[i] - r1[i] - r2[i]) / (2.0 * sqrt(3.0));
			gamma[i] = -(2 * r4[i] - r1[i] - r2[i]) / (2.0 * sqrt(3.0));
			r21[i]   = r2[i] - r1[i];
			r31[i]   = r3[i] - r1[i];
			r32[i]   = r3[i] - r2[i];
			r41[i]   = r4[i] - r1[i];
	        r42[i]   = r4[i] - r2[i];
	}
	// Get the coefficients
	coefficients_const ( r1,  r2,  r3,  r4,  ko,  coef_const,  coefm_const );
	coefficients_f1_f1 ( r1,  r2,  r3,  r4,  ko,  coef_f1_f1,  coefm_f1_f1 );
	coefficients_f1_f2 ( r1,  r2,  r3,  r4,  ko,  coef_f1_f2,  coefm_f1_f2 );
	coefficients_f1_f3 ( r1,  r2,  r3,  r4,  ko,  coef_f1_f3,  coefm_f1_f3 );
	coefficients_f2_f1 ( r1,  r2,  r3,  r4,  ko,  coef_f2_f1,  coefm_f2_f1 );
	coefficients_f2_f2 ( r1,  r2,  r3,  r4,  ko,  coef_f2_f2,  coefm_f2_f2 );
	coefficients_f2_f3 ( r1,  r2,  r3,  r4,  ko,  coef_f2_f3,  coefm_f2_f3 );
	coefficients_f3_f1 ( r1,  r2,  r3,  r4,  ko,  coef_f3_f1,  coefm_f3_f1 );
	coefficients_f3_f2 ( r1,  r2,  r3,  r4,  ko,  coef_f3_f2,  coefm_f3_f2 );
	coefficients_f3_f3 ( r1,  r2,  r3,  r4,  ko,  coef_f3_f3,  coefm_f3_f3 );
	
     // Initialization of I_
	 for ( int im = 0; im <  6; im++ )
	 {
		 I_const[im] = 0.0;
		 I_f1_f1[im] = 0.0;
		 I_f1_f2[im] = 0.0;
		 I_f1_f3[im] = 0.0;
		 I_f2_f1[im] = 0.0;
		 I_f2_f2[im] = 0.0;
		 I_f2_f3[im] = 0.0;
		 I_f3_f1[im] = 0.0;
		 I_f3_f2[im] = 0.0;
		 I_f3_f3[im] = 0.0;
	 }
	 //
	 for ( int m = 1; m <  7; m++ )
	 {
		 I_theta_const = 0.0;
		 I_theta_f1_f1 = 0.0;
		 I_theta_f1_f2 = 0.0;
		 I_theta_f1_f3 = 0.0;
		 I_theta_f2_f1 = 0.0;
		 I_theta_f2_f2 = 0.0;
		 I_theta_f2_f3 = 0.0;
		 I_theta_f3_f1 = 0.0;
		 I_theta_f3_f2 = 0.0;
		 I_theta_f3_f3 = 0.0;
		 //
		 THETA_limits ( m, &theta_A, &theta_B );
		 //
		 for ( int n_theta = 0 ; n_theta <  N_theta ; n_theta++ )
		 {
			 THETA   = ( (theta_B - theta_A) * 0.5) * z_theta[n_theta] + (theta_B + theta_A) * 0.5;
			 J_theta = (theta_B - theta_A) * 0.5;
			 tPsiA   = sin(THETA) - sqrt(3.0) * cos(THETA);
			 tPsiB   = sin(THETA) + sqrt(3.0) * cos(THETA);
			 //
			 PsiA = atan(tPsiA);
			 PsiB = atan(tPsiB);
			 // Evaluate b0, b1, b2
			 b0 = vector_dot(beta,beta);
			 b1 = 2 * ( vector_dot(alpha,beta) * cos(THETA) + vector_dot(beta,gamma) * sin(THETA));
			 b2 = vector_dot(alpha,alpha)  * pow(cos(THETA),2) + 2 * vector_dot(alpha,gamma) * cos(THETA) * sin(THETA) + vector_dot(gamma,gamma) * pow(sin(THETA),2) ;
		     // Evaluate B1, B2 (B_plus, B_minus)
			 B1 = 2 * ( vector_dot(beta,gamma) * sin(THETA) - vector_dot(alpha,beta) * cos(THETA));
			 B2 = vector_dot(alpha,alpha) * pow(cos(THETA),2) - 2 * vector_dot(alpha,gamma) * cos(THETA) * sin(THETA) + vector_dot(gamma,gamma) * pow(sin(THETA),2) ;
		     //
			 PSI_limits ( m, PsiA, PsiB, &psi_A, &psi_B );	
			 //
			 I_psi_const = 0.0;
			 I_psi_f1_f1 = 0.0;
			 I_psi_f1_f2 = 0.0;
			 I_psi_f1_f3 = 0.0;
			 I_psi_f2_f1 = 0.0;
			 I_psi_f2_f2 = 0.0;
			 I_psi_f2_f3 = 0.0;
			 I_psi_f3_f1 = 0.0;
			 I_psi_f3_f2 = 0.0;
			 I_psi_f3_f3 = 0.0;
			 //
			 for ( int n_psi = 0 ; n_psi <  N_psi ; n_psi++ )
			 {
				 PSI   = ( (psi_B - psi_A) * 0.5) * z_psi[n_psi] + (psi_B + psi_A) * 0.5;
				 J_psi = (psi_B - psi_A) * 0.5;
				 //				 
				 B     = sqrt(b0 * pow(sin(PSI),2) + b1 * cos(PSI) * sin(PSI) + b2 * pow(cos(PSI),2) );
				 Bm    = sqrt(b0 * pow(sin(PSI),2) + B1 * cos(PSI) * sin(PSI) + B2 * pow(cos(PSI),2) );
				 // Pre-processing: Evaluate N[12], Nm[12]
				 X_function_pre ( THETA, PSI, tPsiA, tPsiB, PsiA, PsiB, B, Bm, N, Nm, ko);
				 // Evaluate the integrands
				 X_const = X_function_const (THETA, PSI, B, Bm, coef_const, coefm_const, N, Nm);
				 X_f1_f1 = X_function_f1_f1 (THETA, PSI, B, Bm, coef_f1_f1, coefm_f1_f1, N, Nm);
				 X_f1_f2 = X_function_f1_f2 (THETA, PSI, B, Bm, coef_f1_f2, coefm_f1_f2, N, Nm);
				 X_f1_f3 = X_function_f1_f3 (THETA, PSI, B, Bm, coef_f1_f3, coefm_f1_f3, N, Nm);
				 X_f2_f1 = X_function_f2_f1 (THETA, PSI, B, Bm, coef_f2_f1, coefm_f2_f1, N, Nm);
				 X_f2_f2 = X_function_f2_f2 (THETA, PSI, B, Bm, coef_f2_f2, coefm_f2_f2, N, Nm);
				 X_f2_f3 = X_function_f2_f3 (THETA, PSI, B, Bm, coef_f2_f3, coefm_f2_f3, N, Nm);
				 X_f3_f1 = X_function_f3_f1 (THETA, PSI, B, Bm, coef_f3_f1, coefm_f3_f1, N, Nm);
				 X_f3_f2 = X_function_f3_f2 (THETA, PSI, B, Bm, coef_f3_f2, coefm_f3_f2, N, Nm);
				 X_f3_f3 = X_function_f3_f3 (THETA, PSI, B, Bm, coef_f3_f3, coefm_f3_f3, N, Nm);
				 // Integrate
				 WPSI = w_psi[n_psi];
				 //
				 I_psi_const +=  WPSI * X_const;
				 I_psi_f1_f1 +=  WPSI * X_f1_f1;
				 I_psi_f1_f2 +=  WPSI * X_f1_f2;
				 I_psi_f1_f3 +=  WPSI * X_f1_f3;
				 I_psi_f2_f1 +=  WPSI * X_f2_f1;
				 I_psi_f2_f2 +=  WPSI * X_f2_f2;
				 I_psi_f2_f3 +=  WPSI * X_f2_f3;
				 I_psi_f3_f1 +=  WPSI * X_f3_f1;
				 I_psi_f3_f2 +=  WPSI * X_f3_f2;
				 I_psi_f3_f3 +=  WPSI * X_f3_f3;
			 }//end for ( int n_psi = 0 ; n_psi <  N_psi ; n_psi++ )
			 //
			 I_psi_const *=  J_psi;
			 I_psi_f1_f1 *=  J_psi;
			 I_psi_f1_f2 *=  J_psi;
			 I_psi_f1_f3 *=  J_psi;
			 I_psi_f2_f1 *=  J_psi;
			 I_psi_f2_f2 *=  J_psi;
			 I_psi_f2_f3 *=  J_psi;
			 I_psi_f3_f1 *=  J_psi;
			 I_psi_f3_f2 *=  J_psi;
			 I_psi_f3_f3 *=  J_psi;
			 //
			 WTHETA = w_theta[n_theta];
			 //
			 I_theta_const +=  WTHETA * I_psi_const;
			 I_theta_f1_f1 +=  WTHETA * I_psi_f1_f1;
			 I_theta_f1_f2 +=  WTHETA * I_psi_f1_f2;
			 I_theta_f1_f3 +=  WTHETA * I_psi_f1_f3;
			 I_theta_f2_f1 +=  WTHETA * I_psi_f2_f1;
			 I_theta_f2_f2 +=  WTHETA * I_psi_f2_f2;
			 I_theta_f2_f3 +=  WTHETA * I_psi_f2_f3;
			 I_theta_f3_f1 +=  WTHETA * I_psi_f3_f1;
			 I_theta_f3_f2 +=  WTHETA * I_psi_f3_f2;
			 I_theta_f3_f3 +=  WTHETA * I_psi_f3_f3;
		 } //end for ( int n_theta = 0 ; n_theta <  N_theta ; n_theta++ )
		 //
		 I_const[m-1] = J_theta * I_theta_const;
		 I_f1_f1[m-1] = J_theta * I_theta_f1_f1;
		 I_f1_f2[m-1] = J_theta * I_theta_f1_f2;
		 I_f1_f3[m-1] = J_theta * I_theta_f1_f3;
		 I_f2_f1[m-1] = J_theta * I_theta_f2_f1;
		 I_f2_f2[m-1] = J_theta * I_theta_f2_f2;
		 I_f2_f3[m-1] = J_theta * I_theta_f2_f3;
		 I_f3_f1[m-1] = J_theta * I_theta_f3_f1;
		 I_f3_f2[m-1] = J_theta * I_theta_f3_f2;
		 I_f3_f3[m-1] = J_theta * I_theta_f3_f3;
	 } //end for ( int m = 1; m <  7; m++ )
	 //
	 Iconst = I_const[0] + I_const[1] + I_const[2] + I_const[3] + I_const[4] + I_const[5];
	 If1_f1 = I_f1_f1[0] + I_f1_f1[1] + I_f1_f1[2] + I_f1_f1[3] + I_f1_f1[4] + I_f1_f1[5];
	 If1_f2 = I_f1_f2[0] + I_f1_f2[1] + I_f1_f2[2] + I_f1_f2[3] + I_f1_f2[4] + I_f1_f2[5];
	 If1_f3 = I_f1_f3[0] + I_f1_f3[1] + I_f1_f3[2] + I_f1_f3[3] + I_f1_f3[4] + I_f1_f3[5];
	 If2_f1 = I_f2_f1[0] + I_f2_f1[1] + I_f2_f1[2] + I_f2_f1[3] + I_f2_f1[4] + I_f2_f1[5];
	 If2_f2 = I_f2_f2[0] + I_f2_f2[1] + I_f2_f2[2] + I_f2_f2[3] + I_f2_f2[4] + I_f2_f2[5];
	 If2_f3 = I_f2_f3[0] + I_f2_f3[1] + I_f2_f3[2] + I_f2_f3[3] + I_f2_f3[4] + I_f2_f3[5];
	 If3_f1 = I_f3_f1[0] + I_f3_f1[1] + I_f3_f1[2] + I_f3_f1[3] + I_f3_f1[4] + I_f3_f1[5];
	 If3_f2 = I_f3_f2[0] + I_f3_f2[1] + I_f3_f2[2] + I_f3_f2[3] + I_f3_f2[4] + I_f3_f2[5];
	 If3_f3 = I_f3_f3[0] + I_f3_f3[1] + I_f3_f3[2] + I_f3_f3[3] + I_f3_f3[4] + I_f3_f3[5];
	 //
	 double PRECOMP = 1.0 / 12.0;
	 // FINAL OUTPUT
	 I_DE[0] = PRECOMP * sqrt(vector_dot(r32,r32)) * sqrt(vector_dot(r41,r41)) * (Iunit * ko * If1_f1 + 4.0 / (Iunit * ko) * Iconst);
	 I_DE[1] = PRECOMP * sqrt(vector_dot(r32,r32)) * sqrt(vector_dot(r42,r42)) * (Iunit * ko * If1_f2 + 4.0 / (Iunit * ko) * Iconst);
	 I_DE[2] = PRECOMP * sqrt(vector_dot(r32,r32)) * sqrt(vector_dot(r21,r21)) * (Iunit * ko * If1_f3 + 4.0 / (Iunit * ko) * Iconst);
	 I_DE[3] = PRECOMP * sqrt(vector_dot(r31,r31)) * sqrt(vector_dot(r41,r41)) * (Iunit * ko * If2_f1 + 4.0 / (Iunit * ko) * Iconst);
	 I_DE[4] = PRECOMP * sqrt(vector_dot(r31,r31)) * sqrt(vector_dot(r42,r42)) * (Iunit * ko * If2_f2 + 4.0 / (Iunit * ko) * Iconst);
	 I_DE[5] = PRECOMP * sqrt(vector_dot(r31,r31)) * sqrt(vector_dot(r21,r21)) * (Iunit * ko * If2_f3 + 4.0 / (Iunit * ko) * Iconst);
	 I_DE[6] = PRECOMP * sqrt(vector_dot(r21,r21)) * sqrt(vector_dot(r41,r41)) * (Iunit * ko * If3_f1 + 4.0 / (Iunit * ko) * Iconst);
	 I_DE[7] = PRECOMP * sqrt(vector_dot(r21,r21)) * sqrt(vector_dot(r42,r42)) * (Iunit * ko * If3_f2 + 4.0 / (Iunit * ko) * Iconst);
	 I_DE[8] = PRECOMP * sqrt(vector_dot(r21,r21)) * sqrt(vector_dot(r21,r21)) * (Iunit * ko * If3_f3 + 4.0 / (Iunit * ko) * Iconst);
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_const
// ***********************************************************************

void coefficients_const (const double r1[],const double r2[],const double r3[],const double r4[],const double ko, complex<double> coef[], complex<double> coefm[] )
{                         
	complex<double> c[1], cm[1];
	//
	complex<double> j   = Iunit;
	//
	c[0] = 1.0;
	//
	complex<double> t1 = j * j;
	complex<double> t4 = ko * ko;
	complex<double> t6 = c[0] / t1 / t4;
	//
	coef[0] = -t6;
	coef[1] = -c[0] / j / ko;
	coef[2] = t6;
	//
	cm[0] = 1.0;
	//
	t1 = j * j;
	t4 = ko * ko;
	t6 = cm[0] / t1 / t4;
	//
	coefm[0] = -t6;
	coefm[1] = t6;
	coefm[2] = -cm[0] / j / ko;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f1_f1
// ***********************************************************************

void coefficients_f1_f1 (const double r1[],const double r2[],const double r3[],const double r4[],const double ko, complex<double> coef[], complex<double> coefm[] )
{                         
	double aa[3], bb[3], cc[3], dd[3];
	complex<double> c[10], cm[10];
	//
	for (int i = 0; i < 3; i++)
	{
			aa[i] = r1[i] - r1[i];
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
	}
	complex<double> j   = Iunit;
	//
	complex<double> t1 = pow(bb[1], 0.2e1);
	complex<double> t2 = pow(bb[2], 0.2e1);
	complex<double> t3 = pow(bb[0], 0.2e1);
	complex<double> t4 = -t1 - t2 - t3;
	complex<double> t5 = bb[0] * dd[0];
	complex<double> t6 = sqrt(0.3e1);
	complex<double> t9 = bb[1] * dd[1];
	complex<double> t12 = bb[2] * dd[2];
	complex<double> t16 = t2 * t6 / 0.12e2;
	complex<double> t18 = t1 * t6 / 0.12e2;
	complex<double> t20 = t3 * t6 / 0.12e2;
	complex<double> t21 = t5 * t6 / 0.6e1 + t9 * t6 / 0.6e1 + t12 * t6 / 0.6e1 - t16 - t18 - t20;
	complex<double> t22 = cc[0] * bb[0];
	complex<double> t25 = cc[2] * bb[2];
	complex<double> t28 = cc[1] * bb[1];
	complex<double> t31 = -t22 * t6 / 0.6e1 - t25 * t6 / 0.6e1 + t16 + t20 + t18 - t28 * t6 / 0.6e1;
	complex<double> t47 = cc[0] * dd[0] / 0.3e1 + cc[1] * dd[1] / 0.3e1 - t5 / 0.6e1 - t22 / 0.6e1 - t9 / 0.6e1 - t28 / 0.6e1 + t1 / 0.12e2 - t25 / 0.6e1 - t12 / 0.6e1 + cc[2] * dd[2] / 0.3e1 + t3 / 0.12e2 + t2 / 0.12e2;
	//
	c[0] = t4 / 0.4e1;
	c[1] = t21;
	c[2] = t31;
	c[3] = t4 / 0.4e1;
	c[4] = t4 / 0.4e1;
	c[5] = t21;
	c[6] = t31;
	c[7] = t47;
	c[8] = -t4 / 0.4e1;
	c[9] = -t31;
	//
	t1 = ko * ko;
	t3 = j * j;
	t5 = 0.1e1 / t1 / t3;
	t6 =  t5 *  c[0];
	complex<double> t7 =  t5 *  c[8];
	t12 = 0.1e1 / t1 / ko / t3 / j;
	complex<double> t14 = 2.0 * t12 * c[6];
	t16 = 2.0 * t12 * c[5];
	t18 = 2.0 * t12 * c[4];
	t21 = 0.1e1 / ko / j;
	complex<double> t42 = 2.0 * t12 * c[9];
	complex<double> t43 = t1 * t1;
	complex<double> t45 = t3 * t3;
	t47 = 0.1e1 / t43 / t45;
	complex<double> t49 = 6.0 * t47 * c[2];
	complex<double> t51 = 6.0 * t47 * c[7];
	complex<double> t53 = 2.0 * t12 * c[1];
	complex<double> t55 = 2.0 * t12 * c[3];
	//
	coef[0] = - t6;
	coef[1] =  t7;
	coef[2] = t14;
	coef[3] = t16;
	coef[4] = t18;
	coef[5] = -t21 * c[1];
	coef[6] = -t21 * c[3];
	coef[7] = -t21 * c[9];
	coef[8] = -2.0 * t5 * c[9];
	coef[9] = -t21 * c[2];
	coef[10] = -t21 * c[7];
	coef[11] = -2.0 * t5 * c[1];
	coef[12] = -2.0 * t5 * c[3];
	coef[13] = -6.0 * t12 * c[2];
	coef[14] = -6.0 * t12 * c[7];
	coef[15] = -3.0 * t5 * c[2];
	coef[16] = -3.0 * t5 * c[7];
	coef[17] = t42;
	coef[18] = t49;
	coef[19] = t51;
	coef[20] = t53;
	coef[21] = t55;
	coef[22] = -t7;
	coef[23] = -t21 * c[0];
	coef[24] = -t16;
	coef[25] = -t18;
	coef[26] = -t14;
	coef[27] = t6;
	coef[28] = -2.0 * t5 * c[5];
	coef[29] = -t21 * c[8];
	coef[30] = -t21 * c[4];
	coef[31] = -t49;
	coef[32] = -t55;
	coef[33] = -t21 * c[5];
	coef[34] = -t42;
	coef[35] = -2.0 * t5 * c[6];
	coef[36] = -2.0 * t5 * c[4];
	coef[37] = -t21 * c[6];
	coef[38] = -t53;
	coef[39] = -t51;
	//
	t1 = pow(bb[2], 0.2e1);
	t2 = pow(bb[0], 0.2e1);
	t3 = pow(bb[1], 0.2e1);
	t4 = -t1 - t2 - t3;
	t5 = sqrt(0.3e1);
	t7 = t1 * t5 / 0.12e2;
	complex<double> t8 = bb[1] * dd[1];
	complex<double> t11 = bb[2] * dd[2];
	t14 = bb[0] * dd[0];
	t18 = t2 * t5 / 0.12e2;
	t20 = t3 * t5 / 0.12e2;
	t21 = t7 - t8 * t5 / 0.6e1 - t11 * t5 / 0.6e1 - t14 * t5 / 0.6e1 + t18 + t20;
	t22 = cc[2] * bb[2];
	t25 = cc[0] * bb[0];
	t28 = cc[1] * bb[1];
	t31 = t22 * t5 / 0.6e1 - t7 + t25 * t5 / 0.6e1 - t20 + t28 * t5 / 0.6e1 - t18;
	t47 = -t28 / 0.6e1 - t11 / 0.6e1 - t22 / 0.6e1 + cc[2] * dd[2] / 0.3e1 - t25 / 0.6e1 - t14 / 0.6e1 + cc[0] * dd[0] / 0.3e1 - t8 / 0.6e1 + t2 / 0.12e2 + t3 / 0.12e2 + cc[1] * dd[1] / 0.3e1 + t1 / 0.12e2;
	//
	cm[0] = t4 / 0.4e1;
	cm[1] = t4 / 0.4e1;
	cm[2] = t21;
	cm[3] = t31;
	cm[4] = -t21;
	cm[5] = -t4 / 0.4e1;
	cm[6] = -t31;
	cm[7] = -t31;
	cm[8] = -t4 / 0.4e1;
	cm[9] = t47;
	//
	t1 = ko * ko;
	t4 = j * j;
	t7 = 0.1e1 / t1 / ko / t4 / j;
	t9 = 2.0 * t7 * cm[6];
	t11 = 2.0 * t7 * cm[2];
	complex<double> t13 = 2.0 * t7 * cm[1];
	t14 = t1 * t1;
	t16 = t4 * t4;
	t18 = 0.1e1 / t14 / t16;
	t20 = 6.0 * t18 * cm[9];
	t22 = 6.0 * t18 * cm[3];
	complex<double> t24 = 2.0 * t7 * cm[4];
	complex<double> t26 = 2.0 * t7 * cm[5];
	t28 = 2.0 * t7 * cm[7];
	t31 = 0.1e1 / t1 / t4;
	complex<double> t32 = t31 * cm[8];
	complex<double> t35 = 0.1e1 / ko / j;
	complex<double> t37 = t31 * cm[0];
	//
	coefm[0] = t9;
	coefm[1] = t11;
	coefm[2] = t13;
	coefm[3] = t20;
	coefm[4] = t22;
	coefm[5] = -t24;
	coefm[6] = -t26;
	coefm[7] = -t28;
	coefm[8] = -t32;
	coefm[9] = -t35 * cm[0];
	coefm[10] = -t37;
	coefm[11] = t24;
	coefm[12] = t28;
	coefm[13] = t26;
	coefm[14] = t32;
	coefm[15] = -t35 * cm[2];
	coefm[16] = -t35 * cm[1];
	coefm[17] = -2.0 * t31 * cm[6];
	coefm[18] = -t35 * cm[6];
	coefm[19] = -t35 * cm[9];
	coefm[20] = -t35 * cm[3];
	coefm[21] = -2.0 * t31 * cm[2];
	coefm[22] = -2.0 * t31 * cm[1];
	coefm[23] = -6.0 * t7 * cm[9];
	coefm[24] = -6.0 * t7 * cm[3];
	coefm[25] = -3.0 * t31 * cm[9];
	coefm[26] = -3.0 * t31 * cm[3];
	coefm[27] = t37;
	coefm[28] = -2.0 * t31 * cm[5];
	coefm[29] = -t13;
	coefm[30] = -2.0 * t31 * cm[4];
	coefm[31] = -t9;
	coefm[32] = -t35 * cm[8];
	coefm[33] = -t22;
	coefm[34] = -2.0 * t31 * cm[7];
	coefm[35] = -t20;
	coefm[36] = -t11;
	coefm[37] = -t35 * cm[5];
	coefm[38] = -t35 * cm[7];
	coefm[39] = -t35 * cm[4];

}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f1_f2
// ***********************************************************************

void coefficients_f1_f2 (const double r1[],const double r2[],const double r3[],const double r4[],const double ko, complex<double> coef[], complex<double> coefm[] )
{
	double aa[3], bb[3], cc[3], dd[3];
	complex<double> c[11], cm[11];
	//
	for (int i = 0; i < 3; i++)
	{
			aa[i] = r1[i] - r1[i];
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
	}

	complex<double> j = Iunit;
	//
	complex<double> t1 = pow(bb[1], 0.2e1);
	complex<double> t2 = pow(bb[2], 0.2e1);
	complex<double> t3 = pow(bb[0], 0.2e1);
	complex<double> t4 = t1 + t2 + t3;
	complex<double> t5 = bb[2] * dd[2];
	complex<double> t6 = sqrt(0.3e1);
	complex<double> t9 = bb[0] * dd[0];
	complex<double> t12 = bb[1] * dd[1];
	complex<double> t16 = t1 * t6 / 0.12e2;
	complex<double> t18 = t2 * t6 / 0.12e2;
	complex<double> t20 = t3 * t6 / 0.12e2;
	complex<double> t21 = t5 * t6 / 0.6e1 + t9 * t6 / 0.6e1 + t12 * t6 / 0.6e1 - t16 - t18 - t20;
	complex<double> t22 = cc[0] * bb[0];
	complex<double> t25 = cc[1] * bb[1];
	complex<double> t28 = cc[2] * bb[2];
	complex<double> t31 = -t22 * t6 / 0.6e1 + t20 + t16 - t25 * t6 / 0.6e1 - t28 * t6 / 0.6e1 + t18;
	complex<double> t47 = -t25 / 0.6e1 - t28 / 0.6e1 + t2 / 0.12e2 + cc[1] * dd[1] / 0.3e1 - t9 / 0.6e1 + cc[0] * dd[0] / 0.3e1 - t12 / 0.6e1 + t3 / 0.12e2 - t5 / 0.6e1 - t22 / 0.6e1 + t1 / 0.12e2 + cc[2] * dd[2] / 0.3e1;
	//
	c[0] = t4 / 0.4e1;
	c[1] = t21;
	c[2] = t31;
	c[3] = -t4 / 0.4e1;
	c[4] = t21;
	c[5] = -t4 / 0.4e1;
	c[6] = t4 / 0.2e1;
	c[7] = -t31;
	c[8] = t47;
	c[9] = t4 / 0.4e1;
	c[10] = -t31;
	//
	t1 = ko * ko;
	t3 = j * j;
	t5 = 0.1e1 / t1 / t3;
	t6 = t5 * c[0];
	complex<double> t11 = 0.1e1 / t1 / ko / t3 / j;
	complex<double> t13 = 2.0 * t11 * c[5];
	complex<double> t14 = t5 * c[6];
	t16 = 2.0 * t11 * c[7];
	t18 = 2.0 * t11 * c[4];
	complex<double> t19 = t1 * t1;
	t21 = t3 * t3;
	complex<double> t23 = 0.1e1 / t19 / t21;
	t25 = 6.0 * t23 * c[2];
	complex<double> t27 = 6.0 * t23 * c[8];
	complex<double> t29 = 2.0 * t11 * c[3];
	t31 = 2.0 * t11 * c[1];
	complex<double> t33 = 2.0 * t11 * c[10];
	complex<double> t34 = t5 * c[9];
	complex<double> t37 = 0.1e1 / ko / j;
	//
	coef[0] = -t6;
	coef[1] = t13;
	coef[2] = t14;
	coef[3] = t16;
	coef[4] = t18;
	coef[5] = t25;
	coef[6] = t27;
	coef[7] = t29;
	coef[8] = t31;
	coef[9] = t33;
	coef[10] = -t18;
	coef[11] = -t13;
	coef[12] = -t16;
	coef[13] = -t34;
	coef[14] = -t14;
	coef[15] = -t37 * c[0];
	coef[16] = -t37 * c[3];
	coef[17] = -t37 * c[1];
	coef[18] = -t37 * c[10];
	coef[19] = -2.0 * t5 * c[10];
	coef[20] = -t37 * c[2];
	coef[21] = -t37 * c[8];
	coef[22] = -2.0 * t5 * c[3];
	coef[23] = -2.0 * t5 * c[1];
	coef[24] = -6.0 * t11 * c[2];
	coef[25] = -6.0 * t11 * c[8];
	coef[26] = -3.0 * t5 * c[2];
	coef[27] = -3.0 * t5 * c[8];
	coef[28] = t34;
	coef[29] = t6;
	coef[30] = -2.0 * t5 * c[7];
	coef[31] = -t37 * c[9];
	coef[32] = -t33;
	coef[33] = -t37 * c[6];
	coef[34] = -t27;
	coef[35] = -t25;
	coef[36] = -t37 * c[7];
	coef[37] = -2.0 * t5 * c[4];
	coef[38] = -t37 * c[4];
	coef[39] = -2.0 * t5 * c[5];
	coef[40] = -t31;
	coef[41] = -t29;
	coef[42] = -t37 * c[5];
	//
	t1 = pow(bb[1], 0.2e1);
	t2 = pow(bb[2], 0.2e1);
	t3 = pow(bb[0], 0.2e1);
	t4 = t1 + t2 + t3;
	t5 = sqrt(0.3e1);
	complex<double> t7 = t1 * t5 / 0.12e2;
	complex<double> t8 = bb[2] * dd[2];
	t12 = t3 * t5 / 0.12e2;
	t13 = bb[0] * dd[0];
	complex<double> t17 = t2 * t5 / 0.12e2;
	t18 = bb[1] * dd[1];
	t21 = t7 - t8 * t5 / 0.6e1 + t12 - t13 * t5 / 0.6e1 + t17 - t18 * t5 / 0.6e1;
	t22 = cc[0] * bb[0];
	t25 = cc[1] * bb[1];
	t28 = cc[2] * bb[2];
	t31 = t22 * t5 / 0.6e1 - t12 - t17 + t25 * t5 / 0.6e1 + t28 * t5 / 0.6e1 - t7;
	t47 = t1 / 0.12e2 + cc[1] * dd[1] / 0.3e1 - t13 / 0.6e1 - t22 / 0.6e1 - t18 / 0.6e1 - t8 / 0.6e1 + cc[2] * dd[2] / 0.3e1 - t25 / 0.6e1 + t3 / 0.12e2 + t2 / 0.12e2 - t28 / 0.6e1 + cc[0] * dd[0] / 0.3e1;
	//
	cm[0] = t4 / 0.4e1;
	cm[1] = -t4 / 0.4e1;
	cm[2] = t21;
	cm[3] = t31;
	cm[4] = -t4 / 0.2e1;
	cm[5] = t4 / 0.4e1;
	cm[6] = -t21;
	cm[7] = t4 / 0.4e1;
	cm[8] = -t31;
	cm[9] = t31;
	cm[10] = t47;
	//
	t1 = ko * ko;
	t3 = j * j;
	t5 = 0.1e1 / t1 / t3;
	t6 = t5 * cm[0];
	t11 = 0.1e1 / t1 / ko / t3 / j;
	t13 = 2.0 * t11 * cm[6];
	complex<double> t15 = 2.0 * t11 * cm[7];
	t17 = 2.0 * t11 * cm[1];
	t19 = 2.0 * t11 * cm[2];
	t20 = t5 * cm[5];
	t21 = t5 * cm[4];
	t23 = 2.0 * t11 * cm[9];
	complex<double> t26 = 0.1e1 / ko / j;
	t28 = t1 * t1;
	complex<double> t30 = t3 * t3;
	complex<double> t32 = 0.1e1 / t28 / t30;
	t34 = 6.0 * t32 * cm[10];
	complex<double> t36 = 6.0 * t32 * cm[3];
	complex<double> t38 = 2.0 * t11 * cm[8];
	//
	coefm[0] = t6;
	coefm[1] = -t13;
	coefm[2] = -t15;
	coefm[3] = t17;
	coefm[4] = t19;
	coefm[5] = -t20;
	coefm[6] = -t21;
	coefm[7] = -t23;
	coefm[8] = -t26 * cm[0];
	coefm[9] = t34;
	coefm[10] = t23;
	coefm[11] = t36;
	coefm[12] = t38;
	coefm[13] = t20;
	coefm[14] = t21;
	coefm[15] = -t26 * cm[1];
	coefm[16] = -t26 * cm[2];
	coefm[17] = -2.0 * t5 * cm[8];
	coefm[18] = -t26 * cm[8];
	coefm[19] = -t26 * cm[10];
	coefm[20] = -t26 * cm[3];
	coefm[21] = -2.0 * t5 * cm[2];
	coefm[22] = -6.0 * t11 * cm[10];
	coefm[23] = -6.0 * t11 * cm[3];
	coefm[24] = -3.0 * t5 * cm[10];
	coefm[25] = -3.0 * t5 * cm[3];
	coefm[26] = -2.0 * t5 * cm[1];
	coefm[27] = t15;
	coefm[28] = -t6;
	coefm[29] = t13;
	coefm[30] = -t26 * cm[7];
	coefm[31] = -t26 * cm[4];
	coefm[32] = -2.0 * t5 * cm[7];
	coefm[33] = -t19;
	coefm[34] = -2.0 * t5 * cm[9];
	coefm[35] = -t34;
	coefm[36] = -t17;
	coefm[37] = -t26 * cm[6];
	coefm[38] = -t26 * cm[5];
	coefm[39] = -2.0 * t5 * cm[6];
	coefm[40] = -t26 * cm[9];
	coefm[41] = -t36;
	coefm[42] = -t38;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f1_f3
// ***********************************************************************

void coefficients_f1_f3 (const double r1[],const double r2[],const double r3[],const double r4[],const double ko, complex<double> coef[], complex<double> coefm[] )
{
	double aa[3], bb[3], cc[3], dd[3];
	complex<double> c[11], cm[11];
	//
	for (int i = 0; i < 3; i++)
	{
			aa[i] = r1[i] - r1[i];
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
	}
	//
	complex<double> j = Iunit;
	//
	complex<double> t1 = bb[1] * dd[1];
	complex<double> t3 = bb[2] * dd[2];
	complex<double> t5 = pow(bb[0], 0.2e1);
	complex<double> t7 = bb[0] * dd[0];
	complex<double> t9 = pow(bb[2], 0.2e1);
	complex<double> t11 = pow(bb[1], 0.2e1);
	complex<double> t14 = cc[2] * bb[2];
	complex<double> t15 = sqrt(0.3e1);
	complex<double> t17 = t14 * t15 / 0.6e1;
	complex<double> t18 = cc[0] * bb[0];
	complex<double> t20 = t18 * t15 / 0.6e1;
	complex<double> t22 = t11 * t15 / 0.12e2;
	complex<double> t24 = t9 * t15 / 0.12e2;
	complex<double> t26 = t5 * t15 / 0.12e2;
	complex<double> t27 = cc[1] * bb[1];
	complex<double> t29 = t27 * t15 / 0.6e1;
	complex<double> t30 = -t17 - t20 + t22 + t24 + t26 - t29;
	complex<double> t31 = -t11 - t5 - t9;
	complex<double> t33 = t3 * t15 / 0.6e1;
	complex<double> t35 = t7 * t15 / 0.6e1;
	complex<double> t37 = t1 * t15 / 0.6e1;
	complex<double> t38 = t33 + t35 + t37 - t26 - t24 - t22;
	complex<double> t49 = t20 - t24 - t26 + t33 + t29 - t22 + t37 - t15 * cc[0] * dd[0] / 0.3e1 + t35 + t17 - t15 * cc[2] * dd[2] / 0.3e1 - t15 * cc[1] * dd[1] / 0.3e1;
	complex<double> t65 = -t7 / 0.6e1 + cc[0] * dd[0] / 0.3e1 - t27 / 0.6e1 - t18 / 0.6e1 + t5 / 0.12e2 + cc[2] * dd[2] / 0.3e1 + t9 / 0.12e2 - t1 / 0.6e1 - t14 / 0.6e1 + cc[1] * dd[1] / 0.3e1 + t11 / 0.12e2 - t3 / 0.6e1;
	//
	c[0] = -t1 / 0.2e1 - t3 / 0.2e1 + t5 / 0.4e1 - t7 / 0.2e1 + t9 / 0.4e1 + t11 / 0.4e1;
	c[1] = t30;
	c[2] = t31 / 0.4e1;
	c[3] = t38;
	c[4] = t31 / 0.4e1;
	c[5] = -t7 / 0.2e1 - t3 / 0.2e1 + t5 / 0.2e1 + t9 / 0.2e1 + t11 / 0.2e1 - t1 / 0.2e1;
	c[6] = t49;
	c[7] = t65;
	c[8] = -t31 / 0.4e1;
	c[9] = -t30;
	c[10] = t38;
	//
	t1 = ko * ko;
	t3 = j * j;
	t5 = 0.1e1 / t1 / t3;
	complex<double> t6 =  t5 *  c[8];
	t11 = 0.1e1 / t1 / ko / t3 / j;
	complex<double> t13 = 2.0 * t11 * c[4];
	t15 = 2.0 * t11 * c[6];
	complex<double> t16 = t5 * c[5];
	t17 = t5 * c[0];
	complex<double> t19 = 2.0 * t11 * c[3];
	t20 = t1 * t1;
	t22 = t3 * t3;
	t24 = 0.1e1 / t20 / t22;
	t26 = 6.0 * t24 * c[1];
	complex<double> t28 = 6.0 * t24 * c[7];
	t30 = 2.0 * t11 * c[2];
	complex<double> t32 = 2.0 * t11 * c[10];
	complex<double> t34 = 2.0 * t11 * c[9];
	t37 = 0.1e1 / ko / j;
	//
	coef[0] = t6;
	coef[1] = t13;
	coef[2] = t15;
	coef[3] = t16;
	coef[4] = -t17;
	coef[5] = t19;
	coef[6] = t26;
	coef[7] = t28;
	coef[8] = t30;
	coef[9] = t32;
	coef[10] = t34;
	coef[11] = -t19;
	coef[12] = -t13;
	coef[13] = -t15;
	coef[14] = -t6;
	coef[15] = -t16;
	coef[16] = -t37 * c[0];
	coef[17] = -t37 * c[2];
	coef[18] = -t37 * c[10];
	coef[19] = -t37 * c[9];
	coef[20] = -2.0 * t5 * c[9];
	coef[21] = -t37 * c[1];
	coef[22] = -t37 * c[7];
	coef[23] = -2.0 * t5 * c[2];
	coef[24] = -2.0 * t5 * c[10];
	coef[25] = -6.0 * t11 * c[1];
	coef[26] = -6.0 * t11 * c[7];
	coef[27] = -3.0 * t5 * c[1];
	coef[28] = -3.0 * t5 * c[7];
	coef[29] = t17;
	coef[30] = -t34;
	coef[31] = -t37 * c[6];
	coef[32] = -2.0 * t5 * c[6];
	coef[33] = -t37 * c[8];
	coef[34] = -t37 * c[5];
	coef[35] = -t26;
	coef[36] = -t28;
	coef[37] = -2.0 * t5 * c[3];
	coef[38] = -t37 * c[4];
	coef[39] = -t32;
	coef[40] = -t37 * c[3];
	coef[41] = -t30;
	coef[42] = -2.0 * t5 * c[4];
	//
	t1 = pow(bb[0], 0.2e1);
	t3 = pow(bb[1], 0.2e1);
	t5 = pow(bb[2], 0.2e1);
	t7 = bb[0] * dd[0];
	t9 = bb[1] * dd[1];
	t11 = bb[2] * dd[2];
	t15 = -t5 - t1 - t3;
	t16 = sqrt(0.3e1);
	t18 = t11 * t16 / 0.6e1;
	t20 = t7 * t16 / 0.6e1;
	t22 = t1 * t16 / 0.12e2;
	t24 = t5 * t16 / 0.12e2;
	t26 = t9 * t16 / 0.6e1;
	t28 = t3 * t16 / 0.12e2;
	t29 = -t18 - t20 + t22 + t24 - t26 + t28;
	t30 = cc[2] * bb[2];
	t32 = t30 * t16 / 0.6e1;
	t33 = cc[1] * bb[1];
	t35 = t33 * t16 / 0.6e1;
	complex<double> t36 = cc[0] * bb[0];
	t38 = t36 * t16 / 0.6e1;
	complex<double> t39 = t32 - t22 + t35 + t38 - t24 - t28;
	t49 = -t24 + t32 - t22 - t28 + t35 + t38 + t26 + t18 + t20 - t16 * cc[1] * dd[1] / 0.3e1 - t16 * cc[2] * dd[2] / 0.3e1 - t16 * cc[0] * dd[0] / 0.3e1;
	t65 = cc[1] * dd[1] / 0.3e1 - t36 / 0.6e1 - t7 / 0.6e1 - t33 / 0.6e1 - t11 / 0.6e1 + cc[2] * dd[2] / 0.3e1 - t9 / 0.6e1 + t1 / 0.12e2 + t3 / 0.12e2 - t30 / 0.6e1 + t5 / 0.12e2 + cc[0] * dd[0] / 0.3e1;
	//
	cm[0] = t1 / 0.4e1 + t3 / 0.4e1 + t5 / 0.4e1 - t7 / 0.2e1 - t9 / 0.2e1 - t11 / 0.2e1;
	cm[1] = -t3 / 0.2e1 - t5 / 0.2e1 + t7 / 0.2e1 + t9 / 0.2e1 + t11 / 0.2e1 - t1 / 0.2e1;
	cm[2] = t15 / 0.4e1;
	cm[3] = t29;
	cm[4] = t39;
	cm[5] = t49;
	cm[6] = -t15 / 0.4e1;
	cm[7] = t65;
	cm[8] = -t29;
	cm[9] = -t15 / 0.4e1;
	cm[10] = -t39;
	//
	t1 = ko * ko;
	t3 = j * j;
	t5 = 0.1e1 / t1 / t3;
	t6 =  t5 *  cm[0];
	t9 = 0.1e1 / ko / j;
	complex<double> t25 = 0.1e1 / t1 / ko / t3 / j;
	t35 = 2.0 * t25 * cm[3];
	complex<double> t41 = 2.0 * t25 * cm[10];
	complex<double> t43 = t1 * t1;
	complex<double> t45 = t3 * t3;
	complex<double> t47 = 0.1e1 / t43 / t45;
	t49 = 6.0 * t47 * cm[7];
	complex<double> t53 = 6.0 * t47 * cm[4];
	complex<double> t58 = 2.0 * t25 * cm[2];
	complex<double> t60 = t5 * cm[1];
	complex<double> t61 = t5 * cm[6];
	complex<double> t63 = 2.0 * t25 * cm[8];
	t65 = 2.0 * t25 * cm[5];
	complex<double> t67 = 2.0 * t25 * cm[9];
	//
	coefm[0] = -t6;
	coefm[1] = -t9 * cm[2];
	coefm[2] = -t9 * cm[3];
	coefm[3] = -2.0 * t5 * cm[10];
	coefm[4] = -t9 * cm[10];
	coefm[5] = -t9 * cm[7];
	coefm[6] = -t9 * cm[4];
	coefm[7] = -2.0 * t5 * cm[2];
	coefm[8] = -2.0 * t5 * cm[3];
	coefm[9] = -6.0 * t25 * cm[7];
	coefm[10] = -6.0 * t25 * cm[4];
	coefm[11] = -3.0 * t5 * cm[7];
	coefm[12] = -3.0 * t5 * cm[4];
	coefm[13] = -t35;
	coefm[14] = -t9 * cm[9];
	coefm[15] = -t9 * cm[1];
	coefm[16] = -2.0 * t5 * cm[5];
	coefm[17] = -t41;
	coefm[18] = -t9 * cm[6];
	coefm[19] = -t49;
	coefm[20] = -2.0 * t5 * cm[9];
	coefm[21] = -t53;
	coefm[22] = -2.0 * t5 * cm[8];
	coefm[23] = -t9 * cm[5];
	coefm[24] = -t58;
	coefm[25] = -t9 * cm[8];
	coefm[26] = t60;
	coefm[27] = t61;
	coefm[28] = t63;
	coefm[29] = t6;
	coefm[30] = t65;
	coefm[31] = t53;
	coefm[32] = -t63;
	coefm[33] = -t67;
	coefm[34] = -t65;
	coefm[35] = -t9 * cm[0];
	coefm[36] = -t60;
	coefm[37] = -t61;
	coefm[38] = t35;
	coefm[39] = t58;
	coefm[40] = t41;
	coefm[41] = t49;
	coefm[42] = t67;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f2_f1
// ***********************************************************************

void coefficients_f2_f1 (const double r1[],const double r2[],const double r3[],const double r4[],const double ko, complex<double> coef[], complex<double> coefm[] )
{
	double aa[3], bb[3], cc[3], dd[3];
	complex<double> c[11], cm[11];
	//
	for (int i = 0; i < 3; i++)
	{
			aa[i] = r1[i] - r1[i];
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
	}

	complex<double> j = Iunit;
	//
	complex<double> t1 = pow(bb[1], 0.2e1);
	complex<double> t2 = pow(bb[2], 0.2e1);
	complex<double> t3 = pow(bb[0], 0.2e1);
	complex<double> t4 = t1 + t2 + t3;
	complex<double> t5 = bb[2] * dd[2];
	complex<double> t6 = sqrt(0.3e1);
	complex<double> t9 = bb[1] * dd[1];
	complex<double> t12 = bb[0] * dd[0];
	complex<double> t16 = t1 * t6 / 0.12e2;
	complex<double> t18 = t3 * t6 / 0.12e2;
	complex<double> t20 = t2 * t6 / 0.12e2;
	complex<double> t21 = t5 * t6 / 0.6e1 + t9 * t6 / 0.6e1 + t12 * t6 / 0.6e1 - t16 - t18 - t20;
	complex<double> t22 = cc[0] * bb[0];
	complex<double> t25 = cc[1] * bb[1];
	complex<double> t28 = cc[2] * bb[2];
	complex<double> t31 = -t22 * t6 / 0.6e1 + t20 - t25 * t6 / 0.6e1 + t18 - t28 * t6 / 0.6e1 + t16;
	complex<double> t47 = -t25 / 0.6e1 - t9 / 0.6e1 - t28 / 0.6e1 + cc[1] * dd[1] / 0.3e1 - t5 / 0.6e1 + t2 / 0.12e2 + t1 / 0.12e2 - t22 / 0.6e1 + cc[0] * dd[0] / 0.3e1 + t3 / 0.12e2 - t12 / 0.6e1 + cc[2] * dd[2] / 0.3e1;
	//
	c[0] = t4 / 0.4e1;
	c[1] = t21;
	c[2] = t31;
	c[3] = -t4 / 0.4e1;
	c[4] = -t21;
	c[5] = -t4 / 0.2e1;
	c[6] = t31;
	c[7] = t47;
	c[8] = t4 / 0.4e1;
	c[9] = t4 / 0.4e1;
	c[10] = -t31;
	//
	t1 = ko * ko;
	t3 = j * j;
	t5 = 0.1e1 / t1 / t3;
	t6 = t5 * c[0];
	complex<double> t11 = 0.1e1 / t1 / ko / t3 / j;
	complex<double> t13 = 2.0 * t11 * c[9];
	complex<double> t15 = 2.0 * t11 * c[4];
	t16 = t5 * c[5];
	complex<double> t17 = t5 * c[8];
	complex<double> t19 = 2.0 * t11 * c[6];
	t20 = t1 * t1;
	t22 = t3 * t3;
	complex<double> t24 = 0.1e1 / t20 / t22;
	complex<double> t26 = 6.0 * t24 * c[2];
	t28 = 6.0 * t24 * c[7];
	complex<double> t30 = 2.0 * t11 * c[3];
	complex<double> t32 = 2.0 * t11 * c[1];
	complex<double> t34 = 2.0 * t11 * c[10];
	complex<double> t37 = 0.1e1 / ko / j;
	//
	coef[0] = -t6;
	coef[1] = t13;
	coef[2] = t15;
	coef[3] = t16;
	coef[4] = t17;
	coef[5] = t19;
	coef[6] = t26;
	coef[7] = t28;
	coef[8] = t30;
	coef[9] = t32;
	coef[10] = t34;
	coef[11] = -t15;
	coef[12] = -t13;
	coef[13] = -t19;
	coef[14] = -t17;
	coef[15] = -t16;
	coef[16] = -t37 * c[0];
	coef[17] = -t37 * c[3];
	coef[18] = -t37 * c[1];
	coef[19] = -2.0 * t5 * c[10];
	coef[20] = -t37 * c[10];
	coef[21] = -t37 * c[2];
	coef[22] = -t37 * c[7];
	coef[23] = -2.0 * t5 * c[3];
	coef[24] = -2.0 * t5 * c[1];
	coef[25] = -6.0 * t11 * c[2];
	coef[26] = -6.0 * t11 * c[7];
	coef[27] = -3.0 * t5 * c[2];
	coef[28] = -3.0 * t5 * c[7];
	coef[29] = t6;
	coef[30] = -t32;
	coef[31] = -2.0 * t5 * c[6];
	coef[32] = -t37 * c[8];
	coef[33] = -t37 * c[5];
	coef[34] = -t34;
	coef[35] = -t37 * c[6];
	coef[36] = -2.0 * t5 * c[4];
	coef[37] = -t26;
	coef[38] = -t28;
	coef[39] = -t30;
	coef[40] = -t37 * c[9];
	coef[41] = -t37 * c[4];
	coef[42] = -2.0 * t5 * c[9];
	//
	t1 = pow(bb[1], 0.2e1);
	t2 = pow(bb[2], 0.2e1);
	t3 = pow(bb[0], 0.2e1);
	t4 = t1 + t2 + t3;
	t5 = bb[1] * dd[1];
	t6 = sqrt(0.3e1);
	t9 = bb[0] * dd[0];
	t12 = bb[2] * dd[2];
	t16 = t1 * t6 / 0.12e2;
	t18 = t3 * t6 / 0.12e2;
	t20 = t2 * t6 / 0.12e2;
	t21 = -t5 * t6 / 0.6e1 - t9 * t6 / 0.6e1 - t12 * t6 / 0.6e1 + t16 + t18 + t20;
	t22 = cc[2] * bb[2];
	t25 = cc[0] * bb[0];
	t28 = cc[1] * bb[1];
	t31 = -t16 + t22 * t6 / 0.6e1 + t25 * t6 / 0.6e1 - t20 + t28 * t6 / 0.6e1 - t18;
	t47 = -t5 / 0.6e1 - t28 / 0.6e1 + t2 / 0.12e2 + t1 / 0.12e2 - t9 / 0.6e1 - t12 / 0.6e1 - t22 / 0.6e1 + t3 / 0.12e2 - t25 / 0.6e1 + cc[2] * dd[2] / 0.3e1 + cc[0] * dd[0] / 0.3e1 + cc[1] * dd[1] / 0.3e1;
	//
	cm[0] = t4 / 0.4e1;
	cm[1] = t4 / 0.2e1;
	cm[2] = -t4 / 0.4e1;
	cm[3] = t21;
	cm[4] = t31;
	cm[5] = t4 / 0.4e1;
	cm[6] = -t31;
	cm[7] = t47;
	cm[8] = -t31;
	cm[9] = -t4 / 0.4e1;
	cm[10] = t21;
	//
	t1 = ko * ko;
	t4 = j * j;
	complex<double> t7 = 0.1e1 / t1 / ko / t4 / j;
	t9 = 2.0 * t7 * cm[6];
	t12 = 0.1e1 / t1 / t4;
	t13 = t12 * cm[5];
	complex<double> t14 = t12 * cm[1];
	t17 = 0.1e1 / ko / j;
	t20 = 2.0 * t7 * cm[9];
	t22 = 2.0 * t7 * cm[8];
	complex<double> t23 = t12 * cm[0];
	t25 = 2.0 * t7 * cm[10];
	complex<double> t51 = 2.0 * t7 * cm[2];
	complex<double> t53 = t1 * t1;
	complex<double> t55 = t4 * t4;
	complex<double> t57 = 0.1e1 / t53 / t55;
	complex<double> t59 = 6.0 * t57 * cm[7];
	complex<double> t61 = 6.0 * t57 * cm[4];
	complex<double> t68 = 2.0 * t7 * cm[3];
	//
	coefm[0] = t9;
	coefm[1] = -t13;
	coefm[2] = -t14;
	coefm[3] = -t17 * cm[0];
	coefm[4] = t20;
	coefm[5] = t13;
	coefm[6] = t22;
	coefm[7] = t23;
	coefm[8] = t14;
	coefm[9] = t25;
	coefm[10] = -t17 * cm[2];
	coefm[11] = -t17 * cm[3];
	coefm[12] = -2.0 * t12 * cm[6];
	coefm[13] = -t17 * cm[6];
	coefm[14] = -t17 * cm[7];
	coefm[15] = -t17 * cm[4];
	coefm[16] = -2.0 * t12 * cm[2];
	coefm[17] = -2.0 * t12 * cm[3];
	coefm[18] = -6.0 * t7 * cm[7];
	coefm[19] = -6.0 * t7 * cm[4];
	coefm[20] = -3.0 * t12 * cm[7];
	coefm[21] = -3.0 * t12 * cm[4];
	coefm[22] = -t23;
	coefm[23] = -2.0 * t12 * cm[9];
	coefm[24] = -t17 * cm[1];
	coefm[25] = -2.0 * t12 * cm[8];
	coefm[26] = -t51;
	coefm[27] = -t17 * cm[9];
	coefm[28] = -t59;
	coefm[29] = -t61;
	coefm[30] = -2.0 * t12 * cm[10];
	coefm[31] = -t17 * cm[8];
	coefm[32] = -t9;
	coefm[33] = -t17 * cm[5];
	coefm[34] = -t17 * cm[10];
	coefm[35] = -t68;
	coefm[36] = -t22;
	coefm[37] = t51;
	coefm[38] = t68;
	coefm[39] = -t25;
	coefm[40] = -t20;
	coefm[41] = t59;
	coefm[42] = t61;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f2_f2
// ***********************************************************************

void coefficients_f2_f2 (const double r1[],const double r2[],const double r3[],const double r4[],const double ko, complex<double> coef[], complex<double> coefm[] )
{
	double aa[3], bb[3], cc[3], dd[3];
	complex<double> c[10], cm[10];
	//
	for (int i = 0; i < 3; i++)
	{
			aa[i] = r1[i] - r1[i];
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
	}

	complex<double> j = Iunit;
	//
	complex<double> t1 = pow(bb[1], 0.2e1);
	complex<double> t2 = pow(bb[2], 0.2e1);
	complex<double> t3 = pow(bb[0], 0.2e1);
	complex<double> t4 = -t1 - t2 - t3;
	complex<double> t5 = cc[1] * bb[1];
	complex<double> t6 = sqrt(0.3e1);
	complex<double> t10 = t2 * t6 / 0.12e2;
	complex<double> t12 = t1 * t6 / 0.12e2;
	complex<double> t13 = cc[0] * bb[0];
	complex<double> t17 = t3 * t6 / 0.12e2;
	complex<double> t18 = cc[2] * bb[2];
	complex<double> t21 = -t5 * t6 / 0.6e1 + t10 + t12 - t13 * t6 / 0.6e1 + t17 - t18 * t6 / 0.6e1;
	complex<double> t22 = bb[1] * dd[1];
	complex<double> t25 = bb[0] * dd[0];
	complex<double> t28 = bb[2] * dd[2];
	complex<double> t31 = -t22 * t6 / 0.6e1 + t17 + t10 + t12 - t25 * t6 / 0.6e1 - t28 * t6 / 0.6e1;
	complex<double> t47 = -t22 / 0.6e1 - t18 / 0.6e1 + cc[2] * dd[2] / 0.3e1 + cc[0] * dd[0] / 0.3e1 + cc[1] * dd[1] / 0.3e1 - t25 / 0.6e1 - t28 / 0.6e1 - t13 / 0.6e1 + t2 / 0.12e2 - t5 / 0.6e1 + t3 / 0.12e2 + t1 / 0.12e2;
	//
	c[0] = t4 / 0.4e1;
	c[1] = t21;
	c[2] = t4 / 0.4e1;
	c[3] = -t4 / 0.4e1;
	c[4] = t31;
	c[5] = -t21;
	c[6] = t47;
	c[7] = -t4 / 0.4e1;
	c[8] = -t21;
	c[9] = -t31;
	//
	t1 = ko * ko;
	t4 = j * j;
	complex<double> t7 = 0.1e1 / t1 / ko / t4 / j;
	complex<double> t9 = 2.0 * t7 * c[4];
	complex<double> t11 = 2.0 * t7 * c[5];
	complex<double> t14 = 0.1e1 / t1 / t4;
	complex<double> t15 = t14 * c[0];
	complex<double> t16 = t14 * c[7];
	t18 = 2.0 * t7 * c[3];
	t21 = 0.1e1 / ko / j;
	complex<double> t42 = 2.0 * t7 * c[8];
	complex<double> t43 = t1 * t1;
	complex<double> t45 = t4 * t4;
	t47 = 0.1e1 / t43 / t45;
	complex<double> t49 = 6.0 * t47 * c[1];
	complex<double> t51 = 6.0 * t47 * c[6];
	complex<double> t53 = 2.0 * t7 * c[9];
	complex<double> t55 = 2.0 * t7 * c[2];
	//
	coef[0] = t9;
	coef[1] = t11;
	coef[2] = -t15;
	coef[3] = t16;
	coef[4] = t18;
	coef[5] = -t21 * c[9];
	coef[6] = -t21 * c[2];
	coef[7] = -2.0 * t14 * c[8];
	coef[8] = -t21 * c[8];
	coef[9] = -t21 * c[1];
	coef[10] = -t21 * c[6];
	coef[11] = -2.0 * t14 * c[9];
	coef[12] = -2.0 * t14 * c[2];
	coef[13] = -6.0 * t7 * c[1];
	coef[14] = -6.0 * t7 * c[6];
	coef[15] = -3.0 * t14 * c[1];
	coef[16] = -3.0 * t14 * c[6];
	coef[17] = t42;
	coef[18] = t49;
	coef[19] = t51;
	coef[20] = t53;
	coef[21] = t55;
	coef[22] = -t11;
	coef[23] = -t16;
	coef[24] = -t21 * c[0];
	coef[25] = -t9;
	coef[26] = -t18;
	coef[27] = t15;
	coef[28] = -t55;
	coef[29] = -t49;
	coef[30] = -2.0 * t14 * c[3];
	coef[31] = -2.0 * t14 * c[4];
	coef[32] = -2.0 * t14 * c[5];
	coef[33] = -t21 * c[5];
	coef[34] = -t53;
	coef[35] = -t21 * c[7];
	coef[36] = -t21 * c[3];
	coef[37] = -t42;
	coef[38] = -t21 * c[4];
	coef[39] = -t51;
	//
	t1 = pow(bb[1], 0.2e1);
	t2 = pow(bb[2], 0.2e1);
	t3 = pow(bb[0], 0.2e1);
	t4 = -t1 - t2 - t3;
	t5 = sqrt(0.3e1);
	t7 = t2 * t5 / 0.12e2;
	t9 = t3 * t5 / 0.12e2;
	t10 = bb[0] * dd[0];
	t13 = bb[2] * dd[2];
	t16 = bb[1] * dd[1];
	complex<double> t20 = t1 * t5 / 0.12e2;
	t21 = t7 + t9 - t10 * t5 / 0.6e1 - t13 * t5 / 0.6e1 - t16 * t5 / 0.6e1 + t20;
	t22 = cc[2] * bb[2];
	t25 = cc[0] * bb[0];
	t28 = cc[1] * bb[1];
	t31 = -t20 + t22 * t5 / 0.6e1 + t25 * t5 / 0.6e1 - t9 + t28 * t5 / 0.6e1 - t7;
	t47 = -t10 / 0.6e1 - t22 / 0.6e1 + cc[2] * dd[2] / 0.3e1 + t1 / 0.12e2 - t13 / 0.6e1 - t25 / 0.6e1 - t16 / 0.6e1 - t28 / 0.6e1 + t2 / 0.12e2 + cc[1] * dd[1] / 0.3e1 + t3 / 0.12e2 + cc[0] * dd[0] / 0.3e1;
	//
	cm[0] = t4 / 0.4e1;
	cm[1] = t4 / 0.4e1;
	cm[2] = t21;
	cm[3] = t31;
	cm[4] = t31;
	cm[5] = t4 / 0.4e1;
	cm[6] = t21;
	cm[7] = t47;
	cm[8] = -t31;
	cm[9] = -t4 / 0.4e1;
	//
	t3 = 0.1e1 / ko / j;
	t6 = ko * ko;
	complex<double> t8 = j * j;
	t10 = 0.1e1 / t6 / t8;
	complex<double> t24 = 0.1e1 / t6 / ko / t8 / j;
	complex<double> t33 = t6 * t6;
	complex<double> t35 = t8 * t8;
	complex<double> t37 = 0.1e1 / t33 / t35;
	complex<double> t39 = 6.0 * t37 * cm[7];
	complex<double> t41 = 6.0 * t37 * cm[3];
	t43 = 2.0 * t24 * cm[4];
	complex<double> t44 = t10 * cm[9];
	t47 = 2.0 * t24 * cm[8];
	t49 = 2.0 * t24 * cm[2];
	t51 = 2.0 * t24 * cm[1];
	t53 = 2.0 * t24 * cm[6];
	t55 = 2.0 * t24 * cm[5];
	complex<double> t66 = t10 * cm[0];
	//
	coefm[0] = -t3 * cm[2];
	coefm[1] = -t3 * cm[1];
	coefm[2] = -2.0 * t10 * cm[8];
	coefm[3] = -t3 * cm[8];
	coefm[4] = -t3 * cm[7];
	coefm[5] = -t3 * cm[3];
	coefm[6] = -2.0 * t10 * cm[2];
	coefm[7] = -2.0 * t10 * cm[1];
	coefm[8] = -6.0 * t24 * cm[7];
	coefm[9] = -6.0 * t24 * cm[3];
	coefm[10] = -3.0 * t10 * cm[7];
	coefm[11] = -3.0 * t10 * cm[3];
	coefm[12] = t39;
	coefm[13] = t41;
	coefm[14] = -t43;
	coefm[15] = -t44;
	coefm[16] = -t3 * cm[0];
	coefm[17] = t47;
	coefm[18] = t49;
	coefm[19] = t51;
	coefm[20] = -t53;
	coefm[21] = -t55;
	coefm[22] = t53;
	coefm[23] = t55;
	coefm[24] = -t3 * cm[9];
	coefm[25] = -t3 * cm[5];
	coefm[26] = -2.0 * t10 * cm[4];
	coefm[27] = -t51;
	coefm[28] = -2.0 * t10 * cm[6];
	coefm[29] = -t41;
	coefm[30] = -t3 * cm[4];
	coefm[31] = -t39;
	coefm[32] = -t3 * cm[6];
	coefm[33] = -t47;
	coefm[34] = -2.0 * t10 * cm[5];
	coefm[35] = -t49;
	coefm[36] = t44;
	coefm[37] = t43;
	coefm[38] = -t66;
	coefm[39] = t66;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f2_f3
// ***********************************************************************

void coefficients_f2_f3 (const double r1[],const double r2[],const double r3[],const double r4[],const double ko, complex<double> coef[], complex<double> coefm[] )
{
	double aa[3], bb[3], cc[3], dd[3];
	complex<double> c[11], cm[11];
	//
	for (int i = 0; i < 3; i++)
	{
			aa[i] = r1[i] - r1[i];
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
	}
	//
	complex<double> j = Iunit;
	//
	complex<double> t1 = bb[1] * dd[1];
	complex<double> t3 = bb[2] * dd[2];
	complex<double> t5 = pow(bb[1], 0.2e1);
	complex<double> t7 = bb[0] * dd[0];
	complex<double> t9 = pow(bb[2], 0.2e1);
	complex<double> t11 = pow(bb[0], 0.2e1);
	complex<double> t14 = cc[1] * bb[1];
	complex<double> t15 = sqrt(0.3e1);
	complex<double> t17 = t14 * t15 / 0.6e1;
	complex<double> t18 = cc[0] * bb[0];
	complex<double> t20 = t18 * t15 / 0.6e1;
	complex<double> t21 = cc[2] * bb[2];
	complex<double> t23 = t21 * t15 / 0.6e1;
	complex<double> t25 = t5 * t15 / 0.12e2;
	complex<double> t27 = t9 * t15 / 0.12e2;
	complex<double> t29 = t11 * t15 / 0.12e2;
	complex<double> t30 = -t17 - t20 - t23 + t25 + t27 + t29;
	complex<double> t31 = -t11 - t5 - t9;
	complex<double> t33 = t3 * t15 / 0.6e1;
	complex<double> t35 = t1 * t15 / 0.6e1;
	complex<double> t37 = t7 * t15 / 0.6e1;
	complex<double> t38 = t25 - t33 - t35 + t27 - t37 + t29;
	complex<double> t49 = t20 - t29 + t35 + t33 - t27 - t25 - t15 * cc[0] * dd[0] / 0.3e1 - t15 * cc[2] * dd[2] / 0.3e1 + t23 + t17 + t37 - t15 * cc[1] * dd[1] / 0.3e1;
	complex<double> t65 = -t14 / 0.6e1 - t21 / 0.6e1 - t3 / 0.6e1 + t11 / 0.12e2 + t5 / 0.12e2 - t18 / 0.6e1 + cc[1] * dd[1] / 0.3e1 + cc[0] * dd[0] / 0.3e1 + cc[2] * dd[2] / 0.3e1 - t7 / 0.6e1 - t1 / 0.6e1 + t9 / 0.12e2;
	//
	c[0] = t1 / 0.2e1 + t3 / 0.2e1 - t5 / 0.4e1 + t7 / 0.2e1 - t9 / 0.4e1 - t11 / 0.4e1;
	c[1] = t30;
	c[2] = t31 / 0.4e1;
	c[3] = -t31 / 0.4e1;
	c[4] = t38;
	c[5] = -t7 / 0.2e1 - t3 / 0.2e1 - t1 / 0.2e1;
	c[6] = t49;
	c[7] = t65;
	c[8] = -t31 / 0.4e1;
	c[9] = -t30;
	c[10] = -t38;
	//
	t1 = ko * ko;
	t3 = j * j;
	t5 = 0.1e1 / t1 / t3;
	complex<double> t6 = t5 * c[5];
	t7 =  t5 * c[8];
	complex<double> t8 = t5 * c[0];
	complex<double> t13 = 0.1e1 / t1 / ko / t3 / j;
	t15 = 2.0 * t13 * c[6];
	t17 = 2.0 * t13 * c[4];
	complex<double> t19 = 2.0 * t13 * c[3];
	complex<double> t22 = 0.1e1 / ko / j;
	complex<double> t24 = t1 * t1;
	complex<double> t26 = t3 * t3;
	complex<double> t28 = 0.1e1 / t24 / t26;
	t30 = 6.0 * t28 * c[1];
	complex<double> t32 = 6.0 * t28 * c[7];
	complex<double> t34 = 2.0 * t13 * c[2];
	complex<double> t36 = 2.0 * t13 * c[10];
	t38 = 2.0 * t13 * c[9];
	//
	coef[0] = t6;
	coef[1] = t7;
	coef[2] = -t8;
	coef[3] = t15;
	coef[4] = t17;
	coef[5] = t19;
	coef[6] = -t22 * c[0];
	coef[7] = t30;
	coef[8] = t32;
	coef[9] = t34;
	coef[10] = t36;
	coef[11] = t38;
	coef[12] = -t17;
	coef[13] = -t19;
	coef[14] = -t15;
	coef[15] = -t7;
	coef[16] = -t6;
	coef[17] = -t22 * c[10];
	coef[18] = -t22 * c[1];
	coef[19] = -2.0 * t5 * c[9];
	coef[20] = -t22 * c[2];
	coef[21] = -t22 * c[9];
	coef[22] = -t22 * c[7];
	coef[23] = -2.0 * t5 * c[2];
	coef[24] = -2.0 * t5 * c[10];
	coef[25] = -6.0 * t13 * c[1];
	coef[26] = -6.0 * t13 * c[7];
	coef[27] = -3.0 * t5 * c[1];
	coef[28] = -3.0 * t5 * c[7];
	coef[29] = t8;
	coef[30] = -t38;
	coef[31] = -t22 * c[8];
	coef[32] = -t32;
	coef[33] = -2.0 * t5 * c[6];
	coef[34] = -t30;
	coef[35] = -t22 * c[5];
	coef[36] = -t22 * c[6];
	coef[37] = -2.0 * t5 * c[4];
	coef[38] = -t36;
	coef[39] = -t34;
	coef[40] = -t22 * c[4];
	coef[41] = -t22 * c[3];
	coef[42] = -2.0 * t5 * c[3];
	//
	t1 = pow(bb[1], 0.2e1);
	t3 = pow(bb[0], 0.2e1);
	t5 = bb[0] * dd[0];
	t7 = bb[1] * dd[1];
	t9 = bb[2] * dd[2];
	t11 = pow(bb[2], 0.2e1);
	t15 = -t11 - t1 - t3;
	complex<double> t16 = cc[0] * bb[0];
	t17 = sqrt(0.3e1);
	t19 = t16 * t17 / 0.6e1;
	t20 = cc[1] * bb[1];
	t22 = t20 * t17 / 0.6e1;
	t24 = t11 * t17 / 0.12e2;
	t25 = cc[2] * bb[2];
	t27 = t25 * t17 / 0.6e1;
	t29 = t3 * t17 / 0.12e2;
	t31 = t1 * t17 / 0.12e2;
	t32 = t19 + t22 - t24 + t27 - t29 - t31;
	t34 = t7 * t17 / 0.6e1;
	t36 = t5 * t17 / 0.6e1;
	t38 = t9 * t17 / 0.6e1;
	complex<double> t39 = -t34 - t36 + t29 + t24 - t38 + t31;
	t49 = -t29 - t24 + t36 - t17 * cc[1] * dd[1] / 0.3e1 - t31 + t34 + t38 + t27 + t22 + t19 - t17 * cc[2] * dd[2] / 0.3e1 - t17 * cc[0] * dd[0] / 0.3e1;
	t65 = -t5 / 0.6e1 - t9 / 0.6e1 + t3 / 0.12e2 + cc[2] * dd[2] / 0.3e1 + t1 / 0.12e2 - t20 / 0.6e1 + cc[1] * dd[1] / 0.3e1 + cc[0] * dd[0] / 0.3e1 - t7 / 0.6e1 + t11 / 0.12e2 - t25 / 0.6e1 - t16 / 0.6e1;
	//
	cm[0] = -t1 / 0.4e1 - t3 / 0.4e1 + t5 / 0.2e1 + t7 / 0.2e1 + t9 / 0.2e1 - t11 / 0.4e1;
	cm[1] = t7 / 0.2e1 + t5 / 0.2e1 + t9 / 0.2e1;
	cm[2] = t15 / 0.4e1;
	cm[3] = t32;
	cm[4] = t39;
	cm[5] = t15 / 0.4e1;
	cm[6] = t39;
	cm[7] = t49;
	cm[8] = -t15 / 0.4e1;
	cm[9] = -t32;
	cm[10] = t65;
	//
	t3 = 0.1e1 / ko / j;
	t7 = ko * ko;
	t9 = j * j;
	t11 = 0.1e1 / t7 / t9;
	t24 = 0.1e1 / t7 / ko / t9 / j;
	t34 = 2.0 * t24 * cm[6];
	t36 = 2.0 * t24 * cm[5];
	t37 = t11 * cm[1];
	t38 = t11 * cm[8];
	complex<double> t40 = 2.0 * t24 * cm[7];
	complex<double> t42 = 2.0 * t24 * cm[4];
	complex<double> t44 = 2.0 * t24 * cm[2];
	complex<double> t47 = 2.0 * t24 * cm[9];
	complex<double> t48 = t7 * t7;
	complex<double> t50 = t9 * t9;
	complex<double> t52 = 0.1e1 / t48 / t50;
	complex<double> t54 = 6.0 * t52 * cm[3];
	complex<double> t56 = 6.0 * t52 * cm[10];
	complex<double> t57 = t11 * cm[0];
	//
	coefm[0] = -t3 * cm[2];
	coefm[1] = -t3 * cm[4];
	coefm[2] = -t3 * cm[9];
	coefm[3] = -2.0 * t11 * cm[9];
	coefm[4] = -t3 * cm[10];
	coefm[5] = -t3 * cm[3];
	coefm[6] = -2.0 * t11 * cm[2];
	coefm[7] = -2.0 * t11 * cm[4];
	coefm[8] = -6.0 * t24 * cm[10];
	coefm[9] = -6.0 * t24 * cm[3];
	coefm[10] = -3.0 * t11 * cm[10];
	coefm[11] = -3.0 * t11 * cm[3];
	coefm[12] = t34;
	coefm[13] = t36;
	coefm[14] = t37;
	coefm[15] = t38;
	coefm[16] = t40;
	coefm[17] = t42;
	coefm[18] = t44;
	coefm[19] = -t40;
	coefm[20] = -t3 * cm[0];
	coefm[21] = -t37;
	coefm[22] = -t38;
	coefm[23] = t47;
	coefm[24] = -t36;
	coefm[25] = -t34;
	coefm[26] = t54;
	coefm[27] = t56;
	coefm[28] = -t57;
	coefm[29] = t57;
	coefm[30] = -t44;
	coefm[31] = -t54;
	coefm[32] = -t3 * cm[7];
	coefm[33] = -t56;
	coefm[34] = -t42;
	coefm[35] = -t3 * cm[6];
	coefm[36] = -t3 * cm[1];
	coefm[37] = -t3 * cm[5];
	coefm[38] = -2.0 * t11 * cm[7];
	coefm[39] = -2.0 * t11 * cm[5];
	coefm[40] = -2.0 * t11 * cm[6];
	coefm[41] = -t3 * cm[8];
	coefm[42] = -t47;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f3_f1
// ***********************************************************************

void coefficients_f3_f1 (const double r1[],const double r2[],const double r3[],const double r4[],const double ko, complex<double> coef[], complex<double> coefm[] )
{
	double aa[3], bb[3], cc[3], dd[3];
	complex<double> c[11], cm[11];
	//
	for (int i = 0; i < 3; i++)
	{
			aa[i] = r1[i] - r1[i];
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
	}
	//
	complex<double> j = Iunit;
	//
	complex<double> t1 = cc[1] * bb[1];
	complex<double> t3 = cc[2] * bb[2];
	complex<double> t5 = pow(bb[0], 0.2e1);
	complex<double> t7 = pow(bb[1], 0.2e1);
	complex<double> t9 = pow(bb[2], 0.2e1);
	complex<double> t11 = cc[0] * bb[0];
	complex<double> t13 = t1 / 0.2e1 + t3 / 0.2e1 - t5 / 0.4e1 - t7 / 0.4e1 - t9 / 0.4e1 + t11 / 0.2e1;
	complex<double> t14 = sqrt(0.3e1);
	complex<double> t16 = t3 * t14 / 0.6e1;
	complex<double> t18 = t7 * t14 / 0.12e2;
	complex<double> t20 = t5 * t14 / 0.12e2;
	complex<double> t22 = t9 * t14 / 0.12e2;
	complex<double> t24 = t1 * t14 / 0.6e1;
	complex<double> t26 = t11 * t14 / 0.6e1;
	complex<double> t27 = -t16 + t18 + t20 + t22 - t24 - t26;
	complex<double> t28 = -t7 - t5 - t9;
	complex<double> t31 = bb[0] * dd[0];
	complex<double> t35 = bb[2] * dd[2];
	complex<double> t37 = bb[1] * dd[1];
	complex<double> t48 = t9 / 0.12e2 - t31 / 0.6e1 + cc[1] * dd[1] / 0.3e1 - t35 / 0.6e1 - t37 / 0.6e1 + cc[0] * dd[0] / 0.3e1 + t7 / 0.12e2 + t5 / 0.12e2 - t11 / 0.6e1 - t3 / 0.6e1 + cc[2] * dd[2] / 0.3e1 - t1 / 0.6e1;
	complex<double> t50 = t37 * t14 / 0.6e1;
	complex<double> t55 = t35 * t14 / 0.6e1;
	complex<double> t63 = t31 * t14 / 0.6e1;
	complex<double> t64 = t24 + t16 + t26 + t50 - t14 * cc[0] * dd[0] / 0.3e1 + t55 - t14 * cc[2] * dd[2] / 0.3e1 - t14 * cc[1] * dd[1] / 0.3e1 - t22 - t18 + t63 - t20;
	//
	c[0] = t13;
	c[1] = t27;
	c[2] = t28 / 0.4e1;
	c[3] = t13;
	c[4] = -t11 / 0.2e1 - t1 / 0.2e1 - t3 / 0.2e1;
	c[5] = t27;
	c[6] = t48;
	c[7] = -t28 / 0.4e1;
	c[8] = t64;
	c[9] = -t27;
	c[10] = -t18 - t22 + t50 - t20 + t55 + t63;
	//
	t1 = ko * ko;
	complex<double> t4 = j * j;
	t7 = 0.1e1 / t1 / ko / t4 / j;
	t9 = 2.0 * t7 * c[5];
	t11 = 2.0 * t7 * c[3];
	t14 = 0.1e1 / t1 / t4;
	complex<double> t15 = t14 * c[4];
	t16 = t14 * c[7];
	complex<double> t17 = t14 * c[0];
	complex<double> t19 = 2.0 * t7 * c[8];
	t20 = t1 * t1;
	t22 = t4 * t4;
	t24 = 0.1e1 / t20 / t22;
	t26 = 6.0 * t24 * c[1];
	t28 = 6.0 * t24 * c[6];
	complex<double> t30 = 2.0 * t7 * c[2];
	complex<double> t32 = 2.0 * t7 * c[10];
	complex<double> t34 = 2.0 * t7 * c[9];
	t37 = 0.1e1 / ko / j;
	//
	coef[0] = t9;
	coef[1] = t11;
	coef[2] = t15;
	coef[3] = t16;
	coef[4] = -t17;
	coef[5] = t19;
	coef[6] = t26;
	coef[7] = t28;
	coef[8] = t30;
	coef[9] = t32;
	coef[10] = t34;
	coef[11] = -t19;
	coef[12] = -t11;
	coef[13] = -t9;
	coef[14] = -t16;
	coef[15] = -t15;
	coef[16] = -t37 * c[0];
	coef[17] = -t37 * c[2];
	coef[18] = -t37 * c[10];
	coef[19] = -t37 * c[9];
	coef[20] = -2.0 * t14 * c[9];
	coef[21] = -t37 * c[1];
	coef[22] = -t37 * c[6];
	coef[23] = -2.0 * t14 * c[2];
	coef[24] = -2.0 * t14 * c[10];
	coef[25] = -6.0 * t7 * c[1];
	coef[26] = -6.0 * t7 * c[6];
	coef[27] = -3.0 * t14 * c[1];
	coef[28] = -3.0 * t14 * c[6];
	coef[29] = t17;
	coef[30] = -t28;
	coef[31] = -2.0 * t14 * c[8];
	coef[32] = -2.0 * t14 * c[5];
	coef[33] = -t37 * c[7];
	coef[34] = -t26;
	coef[35] = -2.0 * t14 * c[3];
	coef[36] = -t37 * c[5];
	coef[37] = -t37 * c[3];
	coef[38] = -t32;
	coef[39] = -t37 * c[4];
	coef[40] = -t30;
	coef[41] = -t34;
	coef[42] = -t37 * c[8];
	//
	t1 = pow(bb[0], 0.2e1);
	t3 = pow(bb[1], 0.2e1);
	t5 = cc[1] * bb[1];
	t7 = cc[2] * bb[2];
	t9 = pow(bb[2], 0.2e1);
	t11 = cc[0] * bb[0];
	t13 = -t1 / 0.4e1 - t3 / 0.4e1 + t5 / 0.2e1 + t7 / 0.2e1 - t9 / 0.4e1 + t11 / 0.2e1;
	t15 = -t1 - t9 - t3;
	t16 = sqrt(0.3e1);
	t18 = t1 * t16 / 0.12e2;
	t20 = t9 * t16 / 0.12e2;
	t22 = t5 * t16 / 0.6e1;
	t24 = t7 * t16 / 0.6e1;
	t26 = t11 * t16 / 0.6e1;
	t28 = t3 * t16 / 0.12e2;
	complex<double> t29 = -t18 - t20 + t22 + t24 + t26 - t28;
	t30 = bb[1] * dd[1];
	t32 = t30 * t16 / 0.6e1;
	complex<double> t33 = bb[0] * dd[0];
	t35 = t33 * t16 / 0.6e1;
	complex<double> t36 = bb[2] * dd[2];
	complex<double> t38 = t36 * t16 / 0.6e1;
	complex<double> t49 = t35 + t22 - t18 - t28 - t16 * cc[0] * dd[0] / 0.3e1 - t20 + t38 + t26 + t32 + t24 - t16 * cc[2] * dd[2] / 0.3e1 - t16 * cc[1] * dd[1] / 0.3e1;
	complex<double> t65 = cc[1] * dd[1] / 0.3e1 + cc[0] * dd[0] / 0.3e1 - t33 / 0.6e1 - t30 / 0.6e1 + t3 / 0.12e2 - t5 / 0.6e1 + t1 / 0.12e2 + cc[2] * dd[2] / 0.3e1 - t36 / 0.6e1 + t9 / 0.12e2 - t7 / 0.6e1 - t11 / 0.6e1;
	//
	cm[0] = t13;
	cm[1] = t11 / 0.2e1 + t7 / 0.2e1 + t5 / 0.2e1;
	cm[2] = t15 / 0.4e1;
	cm[3] = t29;
	cm[4] = -t32 - t35 - t38 + t20 + t18 + t28;
	cm[5] = -t13;
	cm[6] = t49;
	cm[7] = -t29;
	cm[8] = -t29;
	cm[9] = t65;
	cm[10] = -t15 / 0.4e1;
	//
	t1 = ko * ko;
	t3 = j * j;
	t5 = 0.1e1 / t1 / t3;
	complex<double> t8 = t1 * t1;
	complex<double> t10 = t3 * t3;
	complex<double> t12 = 0.1e1 / t8 / t10;
	t14 = 6.0 * t12 * cm[3];
	t17 = 0.1e1 / ko / j;
	complex<double> t23 = 0.1e1 / t1 / ko / t3 / j;
	complex<double> t25 = 2.0 * t23 * cm[7];
	t29 = 2.0 * t23 * cm[2];
	t36 = 6.0 * t12 * cm[9];
	complex<double> t39 = 2.0 * t23 * cm[4];
	complex<double> t41 = 2.0 * t23 * cm[6];
	complex<double> t42 = t5 * cm[1];
	complex<double> t44 = 2.0 * t23 * cm[8];
	t64 = t5 * cm[0];
	complex<double> t66 = 2.0 * t23 * cm[5];
	complex<double> t67 = t5 * cm[10];
	//
	coefm[0] = -2.0 * t5 * cm[8];
	coefm[1] = -t14;
	coefm[2] = -t17 * cm[5];
	coefm[3] = -t25;
	coefm[4] = -t17 * cm[1];
	coefm[5] = -t17 * cm[6];
	coefm[6] = -t29;
	coefm[7] = -2.0 * t5 * cm[5];
	coefm[8] = -2.0 * t5 * cm[6];
	coefm[9] = -t17 * cm[8];
	coefm[10] = -t36;
	coefm[11] = -t17 * cm[10];
	coefm[12] = -t39;
	coefm[13] = t41;
	coefm[14] = t42;
	coefm[15] = t44;
	coefm[16] = -t17 * cm[2];
	coefm[17] = -t17 * cm[4];
	coefm[18] = -t17 * cm[7];
	coefm[19] = -2.0 * t5 * cm[7];
	coefm[20] = -t17 * cm[9];
	coefm[21] = -t17 * cm[3];
	coefm[22] = -2.0 * t5 * cm[2];
	coefm[23] = -2.0 * t5 * cm[4];
	coefm[24] = -6.0 * t23 * cm[9];
	coefm[25] = -6.0 * t23 * cm[3];
	coefm[26] = -3.0 * t5 * cm[9];
	coefm[27] = -3.0 * t5 * cm[3];
	coefm[28] = t64;
	coefm[29] = t66;
	coefm[30] = t67;
	coefm[31] = t25;
	coefm[32] = -t17 * cm[0];
	coefm[33] = -t67;
	coefm[34] = -t42;
	coefm[35] = -t44;
	coefm[36] = t29;
	coefm[37] = t39;
	coefm[38] = t36;
	coefm[39] = t14;
	coefm[40] = -t41;
	coefm[41] = -t66;
	coefm[42] = -t64;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f3_f2
// ***********************************************************************

void coefficients_f3_f2 (const double r1[],const double r2[],const double r3[],const double r4[],const double ko, complex<double> coef[], complex<double> coefm[] )
{
	double aa[3], bb[3], cc[3], dd[3];
	complex<double> c[11], cm[11];
	//
	for (int i = 0; i < 3; i++)
	{
			aa[i] = r1[i] - r1[i];
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
	}
	//
	complex<double> j = Iunit;
	//
	complex<double> t1 = cc[2] * bb[2];
	complex<double> t3 = pow(bb[0], 0.2e1);
	complex<double> t5 = pow(bb[1], 0.2e1);
	complex<double> t7 = pow(bb[2], 0.2e1);
	complex<double> t9 = cc[1] * bb[1];
	complex<double> t11 = cc[0] * bb[0];
	complex<double> t13 = -t1 / 0.2e1 + t3 / 0.4e1 + t5 / 0.4e1 + t7 / 0.4e1 - t9 / 0.2e1 - t11 / 0.2e1;
	complex<double> t14 = bb[0] * dd[0];
	complex<double> t15 = sqrt(0.3e1);
	complex<double> t17 = t14 * t15 / 0.6e1;
	complex<double> t19 = t7 * t15 / 0.12e2;
	complex<double> t21 = t3 * t15 / 0.12e2;
	complex<double> t26 = t9 * t15 / 0.6e1;
	complex<double> t28 = t5 * t15 / 0.12e2;
	complex<double> t30 = t1 * t15 / 0.6e1;
	complex<double> t32 = t11 * t15 / 0.6e1;
	complex<double> t33 = bb[1] * dd[1];
	complex<double> t35 = t33 * t15 / 0.6e1;
	complex<double> t36 = bb[2] * dd[2];
	complex<double> t38 = t36 * t15 / 0.6e1;
	complex<double> t45 = t17 - t19 - t21 - t15 * cc[0] * dd[0] / 0.3e1 + t26 - t28 + t30 + t32 + t35 + t38 - t15 * cc[2] * dd[2] / 0.3e1 - t15 * cc[1] * dd[1] / 0.3e1;
	complex<double> t47 = t32 + t26 - t19 + t30 - t21 - t28;
	complex<double> t63 = -t9 / 0.6e1 + t7 / 0.12e2 - t1 / 0.6e1 + cc[2] * dd[2] / 0.3e1 - t33 / 0.6e1 + cc[1] * dd[1] / 0.3e1 - t36 / 0.6e1 - t14 / 0.6e1 - t11 / 0.6e1 + t3 / 0.12e2 + cc[0] * dd[0] / 0.3e1 + t5 / 0.12e2;
	complex<double> t64 = t5 + t7 + t3;
	//
	c[0] = t13;
	c[1] = t45;
	c[2] = -t13;
	c[3] = -t11 / 0.2e1 + t5 / 0.2e1 + t3 / 0.2e1 - t1 / 0.2e1 - t9 / 0.2e1 + t7 / 0.2e1;
	c[4] = t47;
	c[5] = t63;
	c[6] = t64 / 0.4e1;
	c[7] = t47;
	c[8] = -t21 + t35 - t28 + t17 - t19 + t38;
	c[9] = -t47;
	c[10] = -t64 / 0.4e1;
	//
	t1 = ko * ko;
	t3 = j * j;
	t5 = 0.1e1 / t1 / t3;
	complex<double> t6 =  t5 * c[3];
	t11 = 0.1e1 / t1 / ko / t3 / j;
	t13 = 2.0 * t11 * c[1];
	t14 = t5 * c[6];
	complex<double> t16 = 2.0 * t11 * c[4];
	t17 = t5 * c[0];
	t19 = 2.0 * t11 * c[2];
	complex<double> t20 = t1 * t1;
	complex<double> t22 = t3 * t3;
	complex<double> t24 = 0.1e1 / t20 / t22;
	t26 = 6.0 * t24 * c[5];
	t28 = 6.0 * t24 * c[9];
	t30 = 2.0 * t11 * c[10];
	t32 = 2.0 * t11 * c[8];
	complex<double> t34 = 2.0 * t11 * c[7];
	complex<double> t37 = 0.1e1 / ko / j;
	//
	coef[0] = t6;
	coef[1] = t13;
	coef[2] = t14;
	coef[3] = t16;
	coef[4] = -t17;
	coef[5] = t19;
	coef[6] = t26;
	coef[7] = t28;
	coef[8] = t30;
	coef[9] = t32;
	coef[10] = t34;
	coef[11] = -t13;
	coef[12] = -t19;
	coef[13] = -t16;
	coef[14] = -t14;
	coef[15] = -t6;
	coef[16] = -t37 * c[0];
	coef[17] = -t37 * c[10];
	coef[18] = -t37 * c[8];
	coef[19] = -2.0 * t5 * c[7];
	coef[20] = -t37 * c[7];
	coef[21] = -t37 * c[5];
	coef[22] = -t37 * c[9];
	coef[23] = -2.0 * t5 * c[10];
	coef[24] = -2.0 * t5 * c[8];
	coef[25] = -6.0 * t11 * c[5];
	coef[26] = -6.0 * t11 * c[9];
	coef[27] = -3.0 * t5 * c[5];
	coef[28] = -3.0 * t5 * c[9];
	coef[29] = t17;
	coef[30] = -t28;
	coef[31] = -t37 * c[6];
	coef[32] = -t34;
	coef[33] = -t26;
	coef[34] = -t37 * c[4];
	coef[35] = -2.0 * t5 * c[4];
	coef[36] = -t37 * c[3];
	coef[37] = -2.0 * t5 * c[1];
	coef[38] = -2.0 * t5 * c[2];
	coef[39] = -t37 * c[2];
	coef[40] = -t37 * c[1];
	coef[41] = -t32;
	coef[42] = -t30;
	//
	t1 = cc[2] * bb[2];
	t3 = pow(bb[0], 0.2e1);
	t5 = pow(bb[1], 0.2e1);
	t7 = pow(bb[2], 0.2e1);
	t9 = cc[1] * bb[1];
	t11 = cc[0] * bb[0];
	t13 = -t1 / 0.2e1 + t3 / 0.4e1 + t5 / 0.4e1 + t7 / 0.4e1 - t9 / 0.2e1 - t11 / 0.2e1;
	t15 = -t3 - t5 - t7;
	t16 = sqrt(0.3e1);
	complex<double> t18 = t1 * t16 / 0.6e1;
	t20 = t9 * t16 / 0.6e1;
	t22 = t5 * t16 / 0.12e2;
	t24 = t11 * t16 / 0.6e1;
	t26 = t3 * t16 / 0.12e2;
	t28 = t7 * t16 / 0.12e2;
	complex<double> t29 = t18 + t20 - t22 + t24 - t26 - t28;
	t30 = bb[0] * dd[0];
	t32 = t30 * t16 / 0.6e1;
	t33 = bb[1] * dd[1];
	t35 = t33 * t16 / 0.6e1;
	t36 = bb[2] * dd[2];
	t38 = t36 * t16 / 0.6e1;
	complex<double> t49 = t38 - t16 * cc[2] * dd[2] / 0.3e1 - t16 * cc[1] * dd[1] / 0.3e1 + t20 - t16 * cc[0] * dd[0] / 0.3e1 - t22 - t26 + t24 - t28 + t18 + t35 + t32;
	complex<double> t65 = -t33 / 0.6e1 - t11 / 0.6e1 + cc[0] * dd[0] / 0.3e1 - t36 / 0.6e1 + t3 / 0.12e2 - t30 / 0.6e1 + cc[1] * dd[1] / 0.3e1 + t5 / 0.12e2 + t7 / 0.12e2 - t1 / 0.6e1 + cc[2] * dd[2] / 0.3e1 - t9 / 0.6e1;
	//
	cm[0] = t13;
	cm[1] = -t7 / 0.2e1 - t5 / 0.2e1 + t11 / 0.2e1 + t9 / 0.2e1 - t3 / 0.2e1 + t1 / 0.2e1;
	cm[2] = t15 / 0.4e1;
	cm[3] = t29;
	cm[4] = -t32 + t26 - t35 - t38 + t28 + t22;
	cm[5] = t29;
	cm[6] = t49;
	cm[7] = t13;
	cm[8] = -t29;
	cm[9] = -t15 / 0.4e1;
	cm[10] = t65;
	//
	t1 = ko * ko;
	complex<double> t4 = j * j;
	t7 = 0.1e1 / t1 / ko / t4 / j;
	t9 = 2.0 * t7 * cm[7];
	complex<double> t12 = 0.1e1 / t1 / t4;
	t13 = t12 * cm[0];
	t15 = 2.0 * t7 * cm[8];
	t18 = 0.1e1 / ko / j;
	t20 = t1 * t1;
	t22 = t4 * t4;
	t24 = 0.1e1 / t20 / t22;
	t26 = 6.0 * t24 * cm[3];
	t29 = 6.0 * t24 * cm[10];
	t38 = 2.0 * t7 * cm[4];
	complex<double> t42 = 2.0 * t7 * cm[2];
	complex<double> t43 = t12 * cm[1];
	t45 = 2.0 * t7 * cm[5];
	t47 = 2.0 * t7 * cm[6];
	complex<double> t48 = t12 * cm[9];
	//
	coefm[0] = t9;
	coefm[1] = t13;
	coefm[2] = -t15;
	coefm[3] = -t18 * cm[5];
	coefm[4] = -t26;
	coefm[5] = -t18 * cm[1];
	coefm[6] = -t29;
	coefm[7] = -t18 * cm[7];
	coefm[8] = -2.0 * t12 * cm[6];
	coefm[9] = -2.0 * t12 * cm[7];
	coefm[10] = -2.0 * t12 * cm[5];
	coefm[11] = -t38;
	coefm[12] = -t18 * cm[9];
	coefm[13] = -t18 * cm[6];
	coefm[14] = -t42;
	coefm[15] = t43;
	coefm[16] = t45;
	coefm[17] = t47;
	coefm[18] = t48;
	coefm[19] = -t13;
	coefm[20] = -t18 * cm[2];
	coefm[21] = -t18 * cm[4];
	coefm[22] = -2.0 * t12 * cm[8];
	coefm[23] = -t18 * cm[8];
	coefm[24] = -t18 * cm[10];
	coefm[25] = -t18 * cm[3];
	coefm[26] = -2.0 * t12 * cm[2];
	coefm[27] = -2.0 * t12 * cm[4];
	coefm[28] = -6.0 * t7 * cm[10];
	coefm[29] = -6.0 * t7 * cm[3];
	coefm[30] = -3.0 * t12 * cm[10];
	coefm[31] = -3.0 * t12 * cm[3];
	coefm[32] = -t47;
	coefm[33] = -t9;
	coefm[34] = t42;
	coefm[35] = t38;
	coefm[36] = t15;
	coefm[37] = t29;
	coefm[38] = t26;
	coefm[39] = -t48;
	coefm[40] = -t43;
	coefm[41] = -t18 * cm[0];
	coefm[42] = -t45;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f3_f3
// ***********************************************************************

void coefficients_f3_f3 (const double r1[],const double r2[],const double r3[],const double r4[],const double ko, complex<double> coef[], complex<double> coefm[] )
{
	double aa[3], bb[3], cc[3], dd[3];
	complex<double> c[11], cm[11];
	//
	for (int i = 0; i < 3; i++)
	{
			aa[i] = r1[i] - r1[i];
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
	}
	//
	complex<double> j = Iunit;
	//
	complex<double> t1 = bb[0] * dd[0];
	complex<double> t3 = cc[0] * dd[0];
	complex<double> t4 = cc[1] * dd[1];
	complex<double> t5 = cc[2] * dd[2];
	complex<double> t6 = pow(bb[1], 0.2e1);
	complex<double> t7 = t6 / 0.4e1;
	complex<double> t8 = pow(bb[2], 0.2e1);
	complex<double> t9 = t8 / 0.4e1;
	complex<double> t10 = pow(bb[0], 0.2e1);
	complex<double> t11 = t10 / 0.4e1;
	complex<double> t12 = cc[1] * bb[1];
	complex<double> t13 = t12 / 0.2e1;
	complex<double> t14 = bb[2] * dd[2];
	complex<double> t16 = cc[0] * bb[0];
	complex<double> t17 = t16 / 0.2e1;
	complex<double> t18 = bb[1] * dd[1];
	complex<double> t20 = cc[2] * bb[2];
	complex<double> t21 = t20 / 0.2e1;
	complex<double> t22 = -t1 / 0.2e1 + t3 + t4 + t5 + t7 + t9 + t11 - t13 - t14 / 0.2e1 - t17 - t18 / 0.2e1 - t21;
	complex<double> t23 = sqrt(0.3e1);
	complex<double> t25 = t6 * t23 / 0.12e2;
	complex<double> t27 = t10 * t23 / 0.12e2;
	complex<double> t29 = t14 * t23 / 0.6e1;
	complex<double> t31 = t18 * t23 / 0.6e1;
	complex<double> t33 = t1 * t23 / 0.6e1;
	complex<double> t35 = t8 * t23 / 0.12e2;
	complex<double> t38 = t12 * t23 / 0.6e1;
	complex<double> t40 = t20 * t23 / 0.6e1;
	complex<double> t42 = t16 * t23 / 0.6e1;
	complex<double> t43 = -t38 - t40 + t25 + t27 - t42 + t35;
	complex<double> t44 = -t10 - t6 - t8;
	complex<double> t55 = -t35 + t29 - t25 + t31 + t38 + t33 - t27 + t40 - t23 * cc[1] * dd[1] / 0.3e1 - t23 * cc[0] * dd[0] / 0.3e1 - t23 * cc[2] * dd[2] / 0.3e1 + t42;
	complex<double> t68 = -t20 / 0.6e1 + t3 / 0.3e1 + t4 / 0.3e1 + t5 / 0.3e1 - t16 / 0.6e1 - t18 / 0.6e1 + t10 / 0.12e2 + t6 / 0.12e2 - t1 / 0.6e1 - t12 / 0.6e1 + t8 / 0.12e2 - t14 / 0.6e1;
	//
	c[0] = t22;
	c[1] = -t25 - t27 + t29 + t31 + t33 - t35;
	c[2] = t43;
	c[3] = t44 / 0.4e1;
	c[4] = -t1 / 0.2e1 - t16 / 0.2e1 - t20 / 0.2e1 + t10 / 0.2e1 - t18 / 0.2e1 - t12 / 0.2e1 + t6 / 0.2e1 - t14 / 0.2e1 + t8 / 0.2e1;
	c[5] = t55;
	c[6] = t68;
	c[7] = -t44 / 0.4e1;
	c[8] = t17 + t21 - t7 + t13 - t11 - t9;
	c[9] = t55;
	c[10] = -t43;
	//
	t1 = ko * ko;
	t3 = j * j;
	t5 = 0.1e1 / t1 / t3;
	t6 = t5 * c[0];
	t11 = 0.1e1 / t1 / ko / t3 / j;
	t13 = 2.0 * t11 * c[5];
	t14 = t5 * c[4];
	t16 = 2.0 * t11 * c[8];
	t18 = 2.0 * t11 * c[9];
	complex<double> t19 = t5 * c[7];
	t20 = t1 * t1;
	t22 = t3 * t3;
	complex<double> t24 = 0.1e1 / t20 / t22;
	complex<double> t26 = 6.0 * t24 * c[2];
	complex<double> t28 = 6.0 * t24 * c[6];
	complex<double> t30 = 2.0 * t11 * c[3];
	t33 = 0.1e1 / ko / j;
	complex<double> t59 = 2.0 * t11 * c[10];
	complex<double> t65 = 2.0 * t11 * c[1];
	//
	coef[0] = -t6;
	coef[1] = t13;
	coef[2] = t14;
	coef[3] = t16;
	coef[4] = t18;
	coef[5] = t19;
	coef[6] = t26;
	coef[7] = t28;
	coef[8] = t30;
	coef[9] = -t33 * c[3];
	coef[10] = -t33 * c[1];
	coef[11] = -2.0 * t5 * c[10];
	coef[12] = -t33 * c[10];
	coef[13] = -t33 * c[2];
	coef[14] = -t33 * c[6];
	coef[15] = -2.0 * t5 * c[3];
	coef[16] = -2.0 * t5 * c[1];
	coef[17] = -6.0 * t11 * c[2];
	coef[18] = -6.0 * t11 * c[6];
	coef[19] = -3.0 * t5 * c[2];
	coef[20] = -3.0 * t5 * c[6];
	coef[21] = t6;
	coef[22] = -2.0 * t5 * c[5];
	coef[23] = -t33 * c[7];
	coef[24] = -t33 * c[4];
	coef[25] = -t33 * c[5];
	coef[26] = -t59;
	coef[27] = -t26;
	coef[28] = -t28;
	coef[29] = -t33 * c[8];
	coef[30] = -2.0 * t5 * c[9];
	coef[31] = -t33 * c[9];
	coef[32] = -t30;
	coef[33] = -t65;
	coef[34] = -2.0 * t5 * c[8];
	coef[35] = t65;
	coef[36] = t59;
	coef[37] = -t18;
	coef[38] = -t16;
	coef[39] = -t13;
	coef[40] = -t19;
	coef[41] = -t14;
	coef[42] = -t33 * c[0];
	//
	t1 = pow(bb[2], 0.2e1);
	complex<double> t2 = t1 / 0.4e1;
	t3 = pow(bb[1], 0.2e1);
	t4 = t3 / 0.4e1;
	t5 = pow(bb[0], 0.2e1);
	t6 = t5 / 0.4e1;
	t7 = cc[1] * bb[1];
	t8 = t7 / 0.2e1;
	t9 = bb[1] * dd[1];
	t11 = cc[0] * bb[0];
	t12 = t11 / 0.2e1;
	t13 = cc[1] * dd[1];
	t14 = cc[2] * dd[2];
	complex<double> t15 = cc[2] * bb[2];
	t16 = t15 / 0.2e1;
	t17 = cc[0] * dd[0];
	t18 = bb[2] * dd[2];
	t20 = bb[0] * dd[0];
	t22 = t2 + t4 + t6 - t8 - t9 / 0.2e1 - t12 + t13 + t14 - t16 + t17 - t18 / 0.2e1 - t20 / 0.2e1;
	t23 = -t3 - t5 - t1;
	t24 = sqrt(0.3e1);
	t26 = t7 * t24 / 0.6e1;
	t28 = t3 * t24 / 0.12e2;
	t30 = t5 * t24 / 0.12e2;
	complex<double> t32 = t15 * t24 / 0.6e1;
	complex<double> t34 = t11 * t24 / 0.6e1;
	complex<double> t36 = t1 * t24 / 0.12e2;
	complex<double> t37 = t26 - t28 - t30 + t32 + t34 - t36;
	complex<double> t39 = t20 * t24 / 0.6e1;
	complex<double> t41 = t9 * t24 / 0.6e1;
	t43 = t18 * t24 / 0.6e1;
	t55 = -t30 - t36 - t24 * cc[2] * dd[2] / 0.3e1 + t32 + t39 - t24 * cc[1] * dd[1] / 0.3e1 - t24 * cc[0] * dd[0] / 0.3e1 - t28 + t43 + t26 + t34 + t41;
	complex<double> t69 = -t9 / 0.6e1 + t14 / 0.3e1 - t15 / 0.6e1 - t7 / 0.6e1 - t20 / 0.6e1 - t11 / 0.6e1 + t3 / 0.12e2 + t13 / 0.3e1 - t18 / 0.6e1 + t5 / 0.12e2 + t1 / 0.12e2 + t17 / 0.3e1;
	//
	cm[0] = t22;
	cm[1] = t23 / 0.4e1;
	cm[2] = t37;
	cm[3] = -t39 + t28 + t36 + t30 - t41 - t43;
	cm[4] = -t5 / 0.2e1 + t11 / 0.2e1 + t9 / 0.2e1 + t18 / 0.2e1 - t3 / 0.2e1 + t7 / 0.2e1 - t1 / 0.2e1 + t20 / 0.2e1 + t15 / 0.2e1;
	cm[5] = t55;
	cm[6] = -t37;
	cm[7] = -t23 / 0.4e1;
	cm[8] = t55;
	cm[9] = -t16 + t2 + t6 - t12 - t8 + t4;
	cm[10] = t69;
	//
	t1 = ko * ko;
	t4 = j * j;
	t7 = 0.1e1 / t1 / ko / t4 / j;
	t9 = 2.0 * t7 * cm[9];
	t12 = 0.1e1 / ko / j;
	t14 = t1 * t1;
	t16 = t4 * t4;
	t18 = 0.1e1 / t14 / t16;
	t20 = 6.0 * t18 * cm[2];
	t22 = 2.0 * t7 * cm[6];
	t24 = 6.0 * t18 * cm[10];
	t27 = 2.0 * t7 * cm[1];
	t30 = 2.0 * t7 * cm[3];
	t34 = 0.1e1 / t1 / t4;
	t43 = 2.0 * t7 * cm[8];
	complex<double> t45 = 2.0 * t7 * cm[5];
	t65 = t34 * cm[4];
	complex<double> t66 = t34 * cm[0];
	complex<double> t67 = t34 * cm[7];
	//
	coefm[0] = t9;
	coefm[1] = -t12 * cm[5];
	coefm[2] = -t20;
	coefm[3] = -t22;
	coefm[4] = -t24;
	coefm[5] = -t12 * cm[7];
	coefm[6] = -t27;
	coefm[7] = -t12 * cm[9];
	coefm[8] = -t30;
	coefm[9] = -t12 * cm[4];
	coefm[10] = -2.0 * t34 * cm[5];
	coefm[11] = -2.0 * t34 * cm[9];
	coefm[12] = -t12 * cm[8];
	coefm[13] = -2.0 * t34 * cm[8];
	coefm[14] = t43;
	coefm[15] = t45;
	coefm[16] = -t12 * cm[1];
	coefm[17] = -t12 * cm[3];
	coefm[18] = -2.0 * t34 * cm[6];
	coefm[19] = -t12 * cm[6];
	coefm[20] = -t12 * cm[10];
	coefm[21] = -t12 * cm[2];
	coefm[22] = -2.0 * t34 * cm[1];
	coefm[23] = -2.0 * t34 * cm[3];
	coefm[24] = -6.0 * t7 * cm[10];
	coefm[25] = -6.0 * t7 * cm[2];
	coefm[26] = -3.0 * t34 * cm[10];
	coefm[27] = -3.0 * t34 * cm[2];
	coefm[28] = t65;
	coefm[29] = -t66;
	coefm[30] = t67;
	coefm[31] = t66;
	coefm[32] = t24;
	coefm[33] = t20;
	coefm[34] = t22;
	coefm[35] = -t43;
	coefm[36] = t30;
	coefm[37] = t27;
	coefm[38] = -t45;
	coefm[39] = -t12 * cm[0];
	coefm[40] = -t65;
	coefm[41] = -t67;
	coefm[42] = -t9;
}

// ***********************************************************************
//			IMPLEMENTATION OF void PSI_limits
// ***********************************************************************

void PSI_limits ( int argument, double PsiA, double PsiB, double *psi_A, double *psi_B )
{
	switch(argument)
	{
		case 1:
			 *psi_A = PsiB;
			 *psi_B = M_PI_2;
			 break;
		case 2:
			 *psi_A = PsiA;
			 *psi_B = M_PI_2;
			 break;
		case 3:
			 *psi_A = 0;
			 *psi_B = PsiA;
			 break;
		case 4:
			 *psi_A = 0;
			 *psi_B = PsiB;
			 break;
		case 5:
			 *psi_A = PsiA;
			 *psi_B = PsiB;
			 break;
		case 6:
			 *psi_A = 0;
			 *psi_B = PsiA;
			 break;
	}
}

// ***********************************************************************
//			IMPLEMENTATION OF void THETA_limits
// ***********************************************************************

void THETA_limits ( int argument, double *theta_A, double *theta_B )
{
	switch(argument)
	{
		case 1:
			 *theta_A = 0.0;
			 *theta_B = M_PI_2 ;
			 break;
		case 2:
			 *theta_A = M_PI_2;
			 *theta_B = M_PI;
			 break;
		case 3:
			 *theta_A = M_PI_2;
			 *theta_B = M_PI;
			 break;
		case 4:
			 *theta_A = 0.0;
			 *theta_B = M_PI / 3;
			 break;
		case 5:
			 *theta_A = M_PI / 3;
			 *theta_B = M_PI_2;
			 break;
		case 6:
			 *theta_A = M_PI / 3;
			 *theta_B = M_PI_2;
			 break;
	}
}

// ***********************************************************************
//			IMPLEMENTATION OF void void X1
// ***********************************************************************

void X1 ( double psi, complex<double> D1, complex<double> D2, double tpsiB, complex<double> a, complex<double> N[] )
{
	complex<double> D3, aD3, expaD3, aD1, expaD1, D2a, H, T11, T21, T31, T41, T12, T22, T32, T42;
	complex<double>	N21, N31, N81, N91, N101, N111, N121, N41, N51, N22, N32, N82, N92, N102, N112, N122, N42, N52;
	//
	D3     = sqrt(3.0) / (cos(psi) * tpsiB);
	aD3    = a * D3;
	expaD3 = exp(-aD3);
	aD1    = a * D1;
	expaD1 = exp(-aD1);
	D2a    = D2 / a;
	//
	H      = 1.0 - D1 * D2;
	//
	T11    = (expaD3 - expaD1) / aD3;
	T21    = (T11 - H * expaD1) / aD3;
	T31    = (2.0 * T21 - pow(H,2) * expaD1) / aD3;
	T41    = (3.0 * T31 - pow(H,3) * expaD1) / aD3;
	//
	T12    = D2a * (1.0 - expaD1);
	T22    = D2a * (1.0 - H * expaD1 - T12);
	T32    = D2a * (1.0 - pow(H,2) * expaD1 - 2.0 * T22);
	T42    = D2a * (1.0 - pow(H,3) * expaD1 - 3.0 * T32);
	//
	N21    = T11;
	N31    = D3 * (T11 + T21);
	N81    = T21;
	N91    = T31;
	N101   = D3 * (T21 + T31);
	N111   = D3 * (T31 + T41);
	N121   = D3 * (N101 + N111);
	N41    = D3 * (N31 + N101);
	N51    = D3 * (N41 + N121);
	//
	N22    = T12;
	N32    = (T12 - T22) / D2;
	N82    = T22;
	N92    = T32;
	N102   = (T22 - T32) / D2;
	N112   = (T32 - T42) / D2;
	N122   = (N102 - N112) / D2;
	N42    = (N32 - N102) / D2;
	N52    = (N42 - N122) / D2;
	// Final Output
	N[0]   = 1.0;
	N[1]   = N21+N22;
	N[2]   = N31+N32;
	N[3]   = N41+N42;
	N[4]   = N51+N52;
	N[5]   = 0.5;
	N[6]   = 1.0 / 3.0;
	N[7]   = N81+N82;
	N[8]   = N91+N92;
	N[9]   = N101+N102;
	N[10]  = N111+N112;
	N[11]  = N121+N122;
}

// ***********************************************************************
//			IMPLEMENTATION OF void void X2
// ***********************************************************************

void X2 ( complex<double> D, complex<double> a, complex<double> N[] )
{
	complex<double> aD, expaD, T1, T2, T3, T4;
	//
	aD    = a * D;
	expaD = exp(-aD);
	//
	T1    = (1.0 - expaD) / aD;
	T2    = (1.0 - T1) / aD;
	T3    = (1.0 - 2.0 * T2) / aD;
	T4    = (1.0 - 3.0 * T3) / aD;
	// Final Output
	N[0]   = 1.0;
	N[1]   = T1;
	N[2]   = D * (T1 - T2);
	N[3]   = 0.0;
	N[4]   = 0.0;
	N[5]   = 0.5;
	N[6]   = 1.0 / 3.0;
	N[7]   = T2;
	N[8]   = T3;
	N[9]   = D * (T2 - T3);
	N[10]  = D * (T3 - T4);
	N[11]  = 0.0;
	//
	N[11]  = D * (N[9] - N[10]);
	N[3]   = D * (N[2] - N[9]);
	N[4]   = D * (N[3] - N[11]);
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> X_function_const
// ***********************************************************************

complex<double> X_function_const (double theta, double Psi, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[])
{                                 
	complex<double> X;
	//
	double t2 = B * B;
	double t4 = 0.1e1 / t2 / B;
	double t7 = cos(Psi);
	double t22 = Bm * Bm;
	double t24 = 0.1e1 / t22 / Bm;
	//
	X = N[0] * t4 * coef[2] * t7 + N[1] * t4 * t7 * coef[0] + N[2] / t2 * t7 * coef[1] + Nm[0] * t24 * coefm[1] * t7 + Nm[1] * t24 * t7 * coefm[0] + Nm[2] / t22 * t7 * coefm[2];
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> X_function_f1_f1
// ***********************************************************************

complex<double> X_function_f1_f1 (double theta, double Psi, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[])
{                                 
	complex<double> X;
	//
	double t2 = B * B;
	double t4 = 0.1e1 / t2 / B;
	double t7 = cos(Psi);
	double t9 = t2 * t2;
	double t10 = 0.1e1 / t9;
	double t13 = t7 * t7;
	double t14 = cos(theta);
	double t15 = t13 * t14;
	double t18 = 0.1e1 / t9 / B;
	double t21 = sin(Psi);
	double t23 = sin(theta);
	double t24 = t13 * t21 * t23;
	double t28 = t7 * t21;
	double t32 = t13 * t23;
	double t36 = t15 * t21;
	double t58 = 0.1e1 / t2;
	complex<double> t173 = N[0] * (t4 * coef[27] * t7 + t10 * coef[4] * t15 + t18 * coef[19] * t24 + t10 * coef[2] * t28 + t10 * coef[3] * t32 + t18 * coef[18] * t36) + N[5] * (t10 * coef[20] * t32 + t10 * coef[21] * t15 + t10 * coef[17] * t28) + N[6] * t4 * t7 * coef[1] + N[10] * t58 * t7 * coef[29] + N[7] * (t10 * coef[32] * t15 + t10 * coef[38] * t32 + t10 * coef[34] * t28) + N[11] * (t58 * coef[5] * t32 + t58 * coef[6] * t15 + t58 * coef[7] * t28) + N[9] * (t4 * coef[11] * t32 + t4 * coef[12] * t15 + t4 * coef[8] * t28) + N[4] * (t58 * coef[9] * t36 + t58 * coef[10] * t24) + N[2] * (t4 * coef[28] * t32 + t4 * coef[36] * t15 + t10 * coef[13] * t36 + t10 * coef[14] * t24 + t58 * coef[23] * t7 + t4 * coef[35] * t28) + N[3] * (t58 * coef[30] * t15 + t58 * coef[33] * t32 + t4 * coef[15] * t36 + t4 * coef[16] * t24 + t58 * coef[37] * t28) + N[1] * (t10 * coef[24] * t32 + t10 * coef[25] * t15 + t18 * coef[31] * t36 + t18 * coef[39] * t24 + t4 * coef[0] * t7 + t10 * coef[26] * t28) + N[8] * t4 * t7 * coef[22];
	double t175 = Bm * Bm;
	double t176 = t175 * t175;
	double t177 = 0.1e1 / t176;
	double t185 = 0.1e1 / t175 / Bm;
	double t190 = 0.1e1 / t176 / Bm;
	double t232 = 0.1e1 / t175;
	complex<double> t335 = Nm[0] * (t177 * coefm[12] * t28 + t177 * coefm[13] * t15 + t185 * coefm[27] * t7 + t190 * coefm[3] * t24 + t190 * coefm[4] * t36 + t177 * coefm[11] * t32) + Nm[5] * (t177 * coefm[1] * t32 + t177 * coefm[2] * t15 + t177 * coefm[0] * t28) + Nm[8] * t185 * t7 * coefm[8] + Nm[7] * (t177 * coefm[29] * t15 + t177 * coefm[36] * t32 + t177 * coefm[31] * t28) + Nm[10] * t232 * t7 * coefm[32] + Nm[6] * t185 * t7 * coefm[14] + Nm[11] * (t232 * coefm[15] * t32 + t232 * coefm[16] * t15 + t232 * coefm[18] * t28) + Nm[9] * (t185 * coefm[21] * t32 + t185 * coefm[22] * t15 + t185 * coefm[17] * t28) + Nm[4] * (t232 * coefm[19] * t24 + t232 * coefm[20] * t36) + Nm[2] * (t185 * coefm[28] * t15 + t185 * coefm[30] * t32 + t177 * coefm[23] * t24 + t177 * coefm[24] * t36 + t232 * coefm[9] * t7 + t185 * coefm[34] * t28) + Nm[3] * (t232 * coefm[37] * t15 + t232 * coefm[39] * t32 + t185 * coefm[25] * t24 + t185 * coefm[26] * t36 + t232 * coefm[38] * t28) + Nm[1] * (t177 * coefm[5] * t32 + t177 * coefm[6] * t15 + t190 * coefm[33] * t36 + t190 * coefm[35] * t24 + t185 * coefm[10] * t7 + t177 * coefm[7] * t28);
	//
	X = t335 + t173;
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> X_function_f1_f2
// ***********************************************************************

complex<double> X_function_f1_f2 (double theta, double Psi, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[])
{
	complex<double> X;
	//
	double t2 = B * B;
	double t3 = t2 * t2;
	double t4 = 0.1e1 / t3;
	double t7 = cos(Psi);
	double t8 = t7 * t7;
	double t9 = sin(theta);
	double t10 = t8 * t9;
	double t14 = cos(theta);
	double t15 = t8 * t14;
	double t18 = 0.1e1 / t3 / B;
	double t21 = sin(Psi);
	double t22 = t15 * t21;
	double t26 = t7 * t21;
	double t29 = 0.1e1 / t2 / B;
	double t36 = t8 * t21 * t9;
	double t53 = 0.1e1 / t2;
	complex<double> t182 = N[0] * (t4 * coef[4] * t10 + t4 * coef[1] * t15 + t18 * coef[5] * t22 + t4 * coef[3] * t26 + t29 * coef[29] * t7 + t18 * coef[6] * t36) + N[2] * (t29 * coef[37] * t10 + t29 * coef[39] * t15 + t4 * coef[24] * t22 + t4 * coef[25] * t36 + t53 * coef[15] * t7 + t29 * coef[30] * t26) + N[5] * (t4 * coef[7] * t15 + t4 * coef[8] * t10 + t29 * coef[2] * t7 + t4 * coef[9] * t26) + N[6] * t29 * t7 * coef[28] + N[11] * (t53 * coef[16] * t15 + t53 * coef[17] * t10 + t53 * coef[18] * t26) + N[4] * (t53 * coef[20] * t22 + t53 * coef[21] * t36) + N[1] * (t4 * coef[10] * t10 + t4 * coef[11] * t15 + t18 * coef[34] * t36 + t18 * coef[35] * t22 + t29 * coef[0] * t7 + t4 * coef[12] * t26) + N[8] * t29 * t7 * coef[13] + N[10] * t53 * t7 * coef[31] + N[7] * (t4 * coef[40] * t10 + t4 * coef[41] * t15 + t29 * coef[14] * t7 + t4 * coef[32] * t26) + N[9] * (t29 * coef[22] * t15 + t29 * coef[23] * t10 + t53 * coef[33] * t7 + t29 * coef[19] * t26) + N[3] * (t53 * coef[38] * t10 + t53 * coef[42] * t15 + t29 * coef[26] * t22 + t29 * coef[27] * t36 + t53 * coef[36] * t26);
	double t184 = Bm * Bm;
	double t186 = 0.1e1 / t184 / Bm;
	double t190 = t184 * t184;
	double t191 = 0.1e1 / t190;
	double t202 = 0.1e1 / t190 / Bm;
	double t227 = 0.1e1 / t184;
	complex<double> t353 = Nm[0] * (t186 * coefm[0] * t7 + t191 * coefm[27] * t15 + t191 * coefm[29] * t10 + t191 * coefm[10] * t26 + t202 * coefm[9] * t36 + t202 * coefm[11] * t22) + Nm[5] * (t191 * coefm[3] * t15 + t191 * coefm[4] * t10 + t186 * coefm[14] * t7 + t191 * coefm[12] * t26) + Nm[11] * (t227 * coefm[15] * t15 + t227 * coefm[16] * t10 + t227 * coefm[18] * t26) + Nm[6] * t186 * t7 * coefm[13] + Nm[2] * (t186 * coefm[32] * t15 + t186 * coefm[39] * t10 + t191 * coefm[22] * t36 + t191 * coefm[23] * t22 + t227 * coefm[8] * t7 + t186 * coefm[34] * t26) + Nm[9] * (t186 * coefm[21] * t10 + t186 * coefm[26] * t15 + t227 * coefm[31] * t7 + t186 * coefm[17] * t26) + Nm[7] * (t191 * coefm[33] * t10 + t191 * coefm[36] * t15 + t186 * coefm[6] * t7 + t191 * coefm[42] * t26) + Nm[10] * t227 * t7 * coefm[38] + Nm[8] * t186 * t7 * coefm[5] + Nm[4] * (t227 * coefm[19] * t36 + t227 * coefm[20] * t22) + Nm[1] * (t191 * coefm[1] * t10 + t191 * coefm[2] * t15 + t202 * coefm[35] * t36 + t202 * coefm[41] * t22 + t186 * coefm[28] * t7 + t191 * coefm[7] * t26) + Nm[3] * (t227 * coefm[30] * t15 + t227 * coefm[37] * t10 + t186 * coefm[24] * t36 + t186 * coefm[25] * t22 + t227 * coefm[40] * t26);
	//
	X = t353 + t182;
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> X_function_f1_f3
// ***********************************************************************

complex<double> X_function_f1_f3 (double theta, double Psi, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[])
{                                 
	complex<double> X;
	//
	double t2 = B * B;
	double t3 = t2 * t2;
	double t4 = 0.1e1 / t3;
	double t7 = cos(Psi);
	double t8 = t7 * t7;
	double t9 = cos(theta);
	double t10 = t8 * t9;
	double t14 = sin(Psi);
	double t15 = t7 * t14;
	double t18 = 0.1e1 / t3 / B;
	double t21 = t10 * t14;
	double t26 = sin(theta);
	double t27 = t8 * t14 * t26;
	double t31 = t8 * t26;
	double t34 = 0.1e1 / t2 / B;
	double t73 = 0.1e1 / t2;
	complex<double> t182 = N[0] * (t4 * coef[1] * t10 + t4 * coef[2] * t15 + t18 * coef[6] * t21 + t18 * coef[7] * t27 + t4 * coef[5] * t31 + t34 * coef[29] * t7) + N[5] * (t4 * coef[8] * t10 + t4 * coef[9] * t31 + t34 * coef[3] * t7 + t4 * coef[10] * t15) + N[6] * t34 * t7 * coef[0] + N[2] * (t34 * coef[37] * t31 + t34 * coef[42] * t10 + t4 * coef[25] * t21 + t4 * coef[26] * t27 + t73 * coef[16] * t7 + t34 * coef[32] * t15) + N[3] * (t73 * coef[38] * t10 + t73 * coef[40] * t31 + t34 * coef[27] * t21 + t34 * coef[28] * t27 + t73 * coef[31] * t15) + N[11] * (t73 * coef[17] * t10 + t73 * coef[18] * t31 + t73 * coef[19] * t15) + N[4] * (t73 * coef[21] * t21 + t73 * coef[22] * t27) + N[1] * (t4 * coef[11] * t31 + t4 * coef[12] * t10 + t18 * coef[35] * t21 + t18 * coef[36] * t27 + t34 * coef[4] * t7 + t4 * coef[13] * t15) + N[8] * t34 * t7 * coef[14] + N[10] * t73 * t7 * coef[33] + N[9] * (t34 * coef[23] * t10 + t34 * coef[24] * t31 + t73 * coef[34] * t7 + t34 * coef[20] * t15) + N[7] * (t4 * coef[39] * t31 + t4 * coef[41] * t10 + t34 * coef[15] * t7 + t4 * coef[30] * t15);
	double t184 = Bm * Bm;
	double t185 = t184 * t184;
	double t187 = 0.1e1 / t185 / Bm;
	double t191 = 0.1e1 / t185;
	double t202 = 0.1e1 / t184 / Bm;
	double t247 = 0.1e1 / t184;
	complex<double> t353 = Nm[0] * (t187 * coefm[31] * t21 + t191 * coefm[28] * t31 + t191 * coefm[30] * t15 + t191 * coefm[42] * t10 + t202 * coefm[29] * t7 + t187 * coefm[41] * t27) + Nm[5] * (t191 * coefm[38] * t31 + t191 * coefm[39] * t10 + t202 * coefm[26] * t7 + t191 * coefm[40] * t15) + Nm[7] * (t191 * coefm[13] * t31 + t191 * coefm[24] * t10 + t202 * coefm[36] * t7 + t191 * coefm[17] * t15) + Nm[8] * t202 * t7 * coefm[37] + Nm[11] * (t247 * coefm[1] * t10 + t247 * coefm[2] * t31 + t247 * coefm[4] * t15) + Nm[6] * t202 * t7 * coefm[27] + Nm[2] * (t202 * coefm[20] * t10 + t202 * coefm[22] * t31 + t191 * coefm[9] * t27 + t191 * coefm[10] * t21 + t247 * coefm[35] * t7 + t202 * coefm[16] * t15) + Nm[3] * (t247 * coefm[14] * t10 + t247 * coefm[25] * t31 + t202 * coefm[11] * t27 + t202 * coefm[12] * t21 + t247 * coefm[23] * t15) + Nm[9] * (t202 * coefm[7] * t10 + t202 * coefm[8] * t31 + t247 * coefm[15] * t7 + t202 * coefm[3] * t15) + Nm[10] * t247 * t7 * coefm[18] + Nm[4] * (t247 * coefm[5] * t27 + t247 * coefm[6] * t21) + Nm[1] * (t191 * coefm[32] * t31 + t191 * coefm[33] * t10 + t187 * coefm[19] * t27 + t187 * coefm[21] * t21 + t202 * coefm[0] * t7 + t191 * coefm[34] * t15);
	//
	X = t353 + t182;
	// Final Output
		return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> X_function_f2_f1
// ***********************************************************************

complex<double> X_function_f2_f1 (double theta, double Psi, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[])
{
	complex<double> X;
	//
	double t2 = B * B;
	double t4 = 0.1e1 / t2 / B;
	double t7 = cos(Psi);
	double t9 = t2 * t2;
	double t10 = 0.1e1 / t9;
	double t13 = t7 * t7;
	double t14 = cos(theta);
	double t15 = t13 * t14;
	double t19 = sin(theta);
	double t20 = t13 * t19;
	double t23 = 0.1e1 / t9 / B;
	double t26 = sin(Psi);
	double t28 = t13 * t26 * t19;
	double t32 = t7 * t26;
	double t36 = t15 * t26;
	double t61 = 0.1e1 / t2;
	complex<double> t182 = N[0] * (t4 * coef[29] * t7 + t10 * coef[1] * t15 + t10 * coef[2] * t20 + t23 * coef[7] * t28 + t10 * coef[5] * t32 + t23 * coef[6] * t36) + N[5] * (t10 * coef[8] * t15 + t10 * coef[9] * t20 + t4 * coef[3] * t7 + t10 * coef[10] * t32) + N[6] * t4 * t7 * coef[4] + N[10] * t61 * t7 * coef[32] + N[9] * (t4 * coef[23] * t15 + t4 * coef[24] * t20 + t61 * coef[33] * t7 + t4 * coef[19] * t32) + N[11] * (t61 * coef[17] * t15 + t61 * coef[18] * t20 + t61 * coef[20] * t32) + N[4] * (t61 * coef[21] * t36 + t61 * coef[22] * t28) + N[1] * (t10 * coef[11] * t20 + t10 * coef[12] * t15 + t23 * coef[37] * t36 + t23 * coef[38] * t28 + t4 * coef[0] * t7 + t10 * coef[13] * t32) + N[8] * t4 * t7 * coef[14] + N[7] * (t10 * coef[30] * t20 + t10 * coef[39] * t15 + t4 * coef[15] * t7 + t10 * coef[34] * t32) + N[3] * (t61 * coef[40] * t15 + t61 * coef[41] * t20 + t4 * coef[27] * t36 + t4 * coef[28] * t28 + t61 * coef[35] * t32) + N[2] * (t4 * coef[36] * t20 + t4 * coef[42] * t15 + t10 * coef[25] * t36 + t10 * coef[26] * t28 + t61 * coef[16] * t7 + t4 * coef[31] * t32);
	double t184 = Bm * Bm;
	double t185 = t184 * t184;
	double t186 = 0.1e1 / t185;
	double t194 = 0.1e1 / t185 / Bm;
	double t202 = 0.1e1 / t184 / Bm;
	double t252 = 0.1e1 / t184;
	complex<double> t353 = Nm[0] * (t186 * coefm[6] * t32 + t186 * coefm[9] * t20 + t194 * coefm[42] * t36 + t186 * coefm[4] * t15 + t202 * coefm[7] * t7 + t194 * coefm[41] * t28) + Nm[5] * (t186 * coefm[37] * t15 + t186 * coefm[38] * t20 + t202 * coefm[8] * t7 + t186 * coefm[0] * t32) + Nm[8] * t202 * t7 * coefm[1] + Nm[7] * (t186 * coefm[26] * t15 + t186 * coefm[35] * t20 + t202 * coefm[2] * t7 + t186 * coefm[32] * t32) + Nm[6] * t202 * t7 * coefm[5] + Nm[11] * (t252 * coefm[10] * t15 + t252 * coefm[11] * t20 + t252 * coefm[13] * t32) + Nm[9] * (t202 * coefm[16] * t15 + t202 * coefm[17] * t20 + t252 * coefm[24] * t7 + t202 * coefm[12] * t32) + Nm[10] * t252 * t7 * coefm[33] + Nm[4] * (t252 * coefm[14] * t28 + t252 * coefm[15] * t36) + Nm[1] * (t186 * coefm[39] * t20 + t186 * coefm[40] * t15 + t194 * coefm[28] * t28 + t194 * coefm[29] * t36 + t202 * coefm[22] * t7 + t186 * coefm[36] * t32) + Nm[3] * (t252 * coefm[27] * t15 + t252 * coefm[34] * t20 + t202 * coefm[20] * t28 + t202 * coefm[21] * t36 + t252 * coefm[31] * t32) + Nm[2] * (t202 * coefm[23] * t15 + t202 * coefm[30] * t20 + t186 * coefm[18] * t28 + t186 * coefm[19] * t36 + t252 * coefm[3] * t7 + t202 * coefm[25] * t32);
	//
	X = t353 + t182;
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> X_function_f2_f2
// ***********************************************************************

complex<double> X_function_f2_f2 (double theta, double Psi, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[])
{
	complex<double> X;
	//
	double t2 = B * B;
	double t4 = 0.1e1 / t2 / B;
	double t7 = cos(Psi);
	double t9 = t2 * t2;
	double t10 = 0.1e1 / t9;
	double t13 = t7 * t7;
	double t14 = sin(theta);
	double t15 = t13 * t14;
	double t19 = sin(Psi);
	double t20 = t7 * t19;
	double t23 = 0.1e1 / t9 / B;
	double t26 = cos(theta);
	double t27 = t13 * t26;
	double t28 = t27 * t19;
	double t36 = t13 * t19 * t14;
	double t58 = 0.1e1 / t2;
	complex<double> t173 = N[0] * (t4 * coef[27] * t7 + t10 * coef[0] * t15 + t10 * coef[1] * t20 + t23 * coef[18] * t28 + t10 * coef[4] * t27 + t23 * coef[19] * t36) + N[5] * (t10 * coef[20] * t15 + t10 * coef[21] * t27 + t10 * coef[17] * t20) + N[6] * t4 * t7 * coef[3] + N[10] * t58 * t7 * coef[35] + N[7] * (t10 * coef[28] * t27 + t10 * coef[34] * t15 + t10 * coef[37] * t20) + N[11] * (t58 * coef[5] * t15 + t58 * coef[6] * t27 + t58 * coef[8] * t20) + N[9] * (t4 * coef[11] * t15 + t4 * coef[12] * t27 + t4 * coef[7] * t20) + N[4] * (t58 * coef[9] * t28 + t58 * coef[10] * t36) + N[2] * (t4 * coef[30] * t27 + t4 * coef[31] * t15 + t10 * coef[13] * t28 + t10 * coef[14] * t36 + t58 * coef[24] * t7 + t4 * coef[32] * t20) + N[3] * (t58 * coef[36] * t27 + t58 * coef[38] * t15 + t4 * coef[15] * t28 + t4 * coef[16] * t36 + t58 * coef[33] * t20) + N[1] * (t10 * coef[25] * t15 + t10 * coef[26] * t27 + t23 * coef[29] * t28 + t23 * coef[39] * t36 + t4 * coef[2] * t7 + t10 * coef[22] * t20) + N[8] * t4 * t7 * coef[23];
	double t175 = Bm * Bm;
	double t176 = t175 * t175;
	double t177 = 0.1e1 / t176;
	double t182 = 0.1e1 / t176 / Bm;
	double t193 = 0.1e1 / t175 / Bm;
	double t220 = 0.1e1 / t175;
	complex<double> t335 = Nm[0] * (t177 * coefm[23] * t27 + t182 * coefm[12] * t36 + t182 * coefm[13] * t28 + t177 * coefm[37] * t20 + t193 * coefm[39] * t7 + t177 * coefm[22] * t15) + Nm[5] * (t177 * coefm[18] * t15 + t177 * coefm[19] * t27 + t177 * coefm[17] * t20) + Nm[8] * t193 * t7 * coefm[15] + Nm[11] * (t220 * coefm[0] * t15 + t220 * coefm[1] * t27 + t220 * coefm[3] * t20) + Nm[9] * (t193 * coefm[6] * t15 + t193 * coefm[7] * t27 + t193 * coefm[2] * t20) + Nm[4] * (t220 * coefm[4] * t36 + t220 * coefm[5] * t28) + Nm[2] * (t193 * coefm[28] * t15 + t193 * coefm[34] * t27 + t177 * coefm[8] * t36 + t177 * coefm[9] * t28 + t220 * coefm[16] * t7 + t193 * coefm[26] * t20) + Nm[3] * (t220 * coefm[25] * t27 + t220 * coefm[32] * t15 + t193 * coefm[10] * t36 + t193 * coefm[11] * t28 + t220 * coefm[30] * t20) + Nm[1] * (t177 * coefm[20] * t15 + t177 * coefm[21] * t27 + t182 * coefm[29] * t28 + t182 * coefm[31] * t36 + t193 * coefm[38] * t7 + t177 * coefm[14] * t20) + Nm[7] * (t177 * coefm[27] * t27 + t177 * coefm[35] * t15 + t177 * coefm[33] * t20) + Nm[10] * t220 * t7 * coefm[24] + Nm[6] * t193 * t7 * coefm[36];
	//
	X = t335 + t173;
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> X_function_f2_f3
// ***********************************************************************

complex<double> X_function_f2_f3 (double theta, double Psi, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[])
{
	complex<double> X;
	//
	double t2 = B * B;
	double t4 = 0.1e1 / t2 / B;
	double t7 = cos(Psi);
	double t9 = t2 * t2;
	double t10 = 0.1e1 / t9;
	double t13 = t7 * t7;
	double t14 = sin(theta);
	double t15 = t13 * t14;
	double t19 = cos(theta);
	double t20 = t13 * t19;
	double t23 = 0.1e1 / t9 / B;
	double t26 = sin(Psi);
	double t27 = t20 * t26;
	double t31 = t7 * t26;
	double t36 = t13 * t26 * t14;
	double t61 = 0.1e1 / t2;
	complex<double> t182 = N[0] * (t4 * coef[29] * t7 + t10 * coef[4] * t15 + t10 * coef[5] * t20 + t23 * coef[7] * t27 + t10 * coef[3] * t31 + t23 * coef[8] * t36) + N[5] * (t10 * coef[9] * t20 + t10 * coef[10] * t15 + t4 * coef[0] * t7 + t10 * coef[11] * t31) + N[6] * t4 * t7 * coef[1] + N[11] * (t61 * coef[17] * t15 + t61 * coef[20] * t20 + t61 * coef[21] * t31) + N[4] * (t61 * coef[18] * t27 + t61 * coef[22] * t36) + N[1] * (t10 * coef[12] * t15 + t10 * coef[13] * t20 + t23 * coef[32] * t36 + t23 * coef[34] * t27 + t4 * coef[2] * t7 + t10 * coef[14] * t31) + N[8] * t4 * t7 * coef[15] + N[10] * t61 * t7 * coef[31] + N[2] * (t4 * coef[37] * t15 + t4 * coef[42] * t20 + t10 * coef[25] * t27 + t10 * coef[26] * t36 + t61 * coef[6] * t7 + t4 * coef[33] * t31) + N[9] * (t4 * coef[23] * t20 + t4 * coef[24] * t15 + t61 * coef[35] * t7 + t4 * coef[19] * t31) + N[3] * (t61 * coef[40] * t15 + t61 * coef[41] * t20 + t4 * coef[27] * t27 + t4 * coef[28] * t36 + t61 * coef[36] * t31) + N[7] * (t10 * coef[38] * t15 + t10 * coef[39] * t20 + t4 * coef[16] * t7 + t10 * coef[30] * t31);
	double t184 = Bm * Bm;
	double t185 = t184 * t184;
	double t186 = 0.1e1 / t185;
	double t191 = 0.1e1 / t185 / Bm;
	double t205 = 0.1e1 / t184 / Bm;
	double t232 = 0.1e1 / t184;
	complex<double> t353 = Nm[0] * (t186 * coefm[13] * t20 + t191 * coefm[26] * t27 + t191 * coefm[27] * t36 + t186 * coefm[16] * t31 + t186 * coefm[12] * t15 + t205 * coefm[29] * t7) + Nm[5] * (t186 * coefm[17] * t15 + t186 * coefm[18] * t20 + t205 * coefm[14] * t7 + t186 * coefm[23] * t31) + Nm[8] * t205 * t7 * coefm[22] + Nm[11] * (t232 * coefm[0] * t20 + t232 * coefm[1] * t15 + t232 * coefm[2] * t31) + Nm[4] * (t232 * coefm[4] * t36 + t232 * coefm[5] * t27) + Nm[1] * (t186 * coefm[24] * t20 + t186 * coefm[25] * t15 + t191 * coefm[31] * t27 + t191 * coefm[33] * t36 + t205 * coefm[28] * t7 + t186 * coefm[19] * t31) + Nm[6] * t205 * t7 * coefm[15] + Nm[7] * (t186 * coefm[30] * t20 + t186 * coefm[34] * t15 + t205 * coefm[21] * t7 + t186 * coefm[42] * t31) + Nm[9] * (t205 * coefm[6] * t20 + t205 * coefm[7] * t15 + t232 * coefm[36] * t7 + t205 * coefm[3] * t31) + Nm[10] * t232 * t7 * coefm[41] + Nm[2] * (t205 * coefm[39] * t20 + t205 * coefm[40] * t15 + t186 * coefm[8] * t36 + t186 * coefm[9] * t27 + t232 * coefm[20] * t7 + t205 * coefm[38] * t31) + Nm[3] * (t232 * coefm[35] * t15 + t232 * coefm[37] * t20 + t205 * coefm[10] * t36 + t205 * coefm[11] * t27 + t232 * coefm[32] * t31);
	//
	X = t353 + t182;
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> X_function_f3_f1
// ***********************************************************************

complex<double> X_function_f3_f1 (double theta, double Psi, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[])
{
	complex<double> X;
	//
	double t2 = B * B;
	double t4 = 0.1e1 / t2 / B;
	double t7 = cos(Psi);
	double t9 = t2 * t2;
	double t10 = 0.1e1 / t9;
	double t13 = sin(Psi);
	double t14 = t7 * t13;
	double t18 = t7 * t7;
	double t19 = cos(theta);
	double t20 = t18 * t19;
	double t23 = 0.1e1 / t9 / B;
	double t27 = sin(theta);
	double t28 = t18 * t13 * t27;
	double t32 = t18 * t27;
	double t36 = t20 * t13;
	double t61 = 0.1e1 / t2;
	complex<double> t182 = N[0] * (t4 * coef[29] * t7 + t10 * coef[0] * t14 + t10 * coef[1] * t20 + t23 * coef[7] * t28 + t10 * coef[5] * t32 + t23 * coef[6] * t36) + N[5] * (t10 * coef[8] * t20 + t10 * coef[9] * t32 + t4 * coef[2] * t7 + t10 * coef[10] * t14) + N[6] * t4 * t7 * coef[3] + N[11] * (t61 * coef[17] * t20 + t61 * coef[18] * t32 + t61 * coef[19] * t14) + N[4] * (t61 * coef[21] * t36 + t61 * coef[22] * t28) + N[1] * (t10 * coef[11] * t32 + t10 * coef[12] * t20 + t23 * coef[30] * t28 + t23 * coef[34] * t36 + t4 * coef[4] * t7 + t10 * coef[13] * t14) + N[8] * t4 * t7 * coef[14] + N[2] * (t4 * coef[31] * t32 + t4 * coef[35] * t20 + t10 * coef[25] * t36 + t10 * coef[26] * t28 + t61 * coef[16] * t7 + t4 * coef[32] * t14) + N[10] * t61 * t7 * coef[33] + N[3] * (t61 * coef[37] * t20 + t61 * coef[42] * t32 + t4 * coef[27] * t36 + t4 * coef[28] * t28 + t61 * coef[36] * t14) + N[7] * (t10 * coef[38] * t32 + t10 * coef[40] * t20 + t4 * coef[15] * t7 + t10 * coef[41] * t14) + N[9] * (t4 * coef[23] * t20 + t4 * coef[24] * t32 + t61 * coef[39] * t7 + t4 * coef[20] * t14);
	double t184 = Bm * Bm;
	double t186 = 0.1e1 / t184 / Bm;
	double t190 = t184 * t184;
	double t191 = 0.1e1 / t190;
	double t199 = 0.1e1 / t190 / Bm;
	double t233 = 0.1e1 / t184;
	complex<double> t353 = Nm[0] * (t186 * coefm[28] * t7 + t191 * coefm[15] * t14 + t191 * coefm[29] * t20 + t199 * coefm[38] * t28 + t199 * coefm[39] * t36 + t191 * coefm[13] * t32) + Nm[5] * (t191 * coefm[36] * t20 + t191 * coefm[37] * t32 + t186 * coefm[14] * t7 + t191 * coefm[31] * t14) + Nm[9] * (t186 * coefm[22] * t20 + t186 * coefm[23] * t32 + t233 * coefm[4] * t7 + t186 * coefm[19] * t14) + Nm[10] * t233 * t7 * coefm[11] + Nm[11] * (t233 * coefm[16] * t20 + t233 * coefm[17] * t32 + t233 * coefm[18] * t14) + Nm[8] * t186 * t7 * coefm[33] + Nm[7] * (t191 * coefm[6] * t20 + t191 * coefm[12] * t32 + t186 * coefm[34] * t7 + t191 * coefm[3] * t14) + Nm[4] * (t233 * coefm[20] * t28 + t233 * coefm[21] * t36) + Nm[1] * (t191 * coefm[40] * t32 + t191 * coefm[41] * t20 + t199 * coefm[1] * t36 + t199 * coefm[10] * t28 + t186 * coefm[42] * t7 + t191 * coefm[35] * t14) + Nm[6] * t186 * t7 * coefm[30] + Nm[2] * (t186 * coefm[7] * t20 + t186 * coefm[8] * t32 + t191 * coefm[24] * t28 + t191 * coefm[25] * t36 + t233 * coefm[32] * t7 + t186 * coefm[0] * t14) + Nm[3] * (t233 * coefm[2] * t20 + t233 * coefm[5] * t32 + t186 * coefm[26] * t28 + t186 * coefm[27] * t36 + t233 * coefm[9] * t14);
	//
	X = t353 + t182;
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> X_function_f3_f2
// ***********************************************************************

complex<double> X_function_f3_f2 (double theta, double Psi, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[])
{
	complex<double> X;
	//
	double t2 = B * B;
	double t4 = 0.1e1 / t2 / B;
	double t7 = cos(Psi);
	double t9 = t2 * t2;
	double t10 = 0.1e1 / t9;
	double t13 = t7 * t7;
	double t14 = sin(theta);
	double t15 = t13 * t14;
	double t18 = 0.1e1 / t9 / B;
	double t21 = cos(theta);
	double t22 = t13 * t21;
	double t23 = sin(Psi);
	double t24 = t22 * t23;
	double t28 = t7 * t23;
	double t36 = t13 * t23 * t14;
	double t82 = 0.1e1 / t2;
	complex<double> t182 = N[0] * (t4 * coef[29] * t7 + t10 * coef[1] * t15 + t18 * coef[7] * t24 + t10 * coef[3] * t28 + t10 * coef[5] * t22 + t18 * coef[6] * t36) + N[5] * (t10 * coef[8] * t22 + t10 * coef[9] * t15 + t4 * coef[0] * t7 + t10 * coef[10] * t28) + N[6] * t4 * t7 * coef[2] + N[1] * (t10 * coef[11] * t15 + t10 * coef[12] * t22 + t18 * coef[30] * t24 + t18 * coef[33] * t36 + t4 * coef[4] * t7 + t10 * coef[13] * t28) + N[11] * (t82 * coef[17] * t22 + t82 * coef[18] * t15 + t82 * coef[20] * t28) + N[4] * (t82 * coef[21] * t36 + t82 * coef[22] * t24) + N[8] * t4 * t7 * coef[14] + N[10] * t82 * t7 * coef[31] + N[7] * (t10 * coef[41] * t15 + t10 * coef[42] * t22 + t4 * coef[15] * t7 + t10 * coef[32] * t28) + N[3] * (t82 * coef[39] * t22 + t82 * coef[40] * t15 + t4 * coef[27] * t36 + t4 * coef[28] * t24 + t82 * coef[34] * t28) + N[2] * (t4 * coef[37] * t15 + t4 * coef[38] * t22 + t10 * coef[25] * t36 + t10 * coef[26] * t24 + t82 * coef[16] * t7 + t4 * coef[35] * t28) + N[9] * (t4 * coef[23] * t22 + t4 * coef[24] * t15 + t82 * coef[36] * t7 + t4 * coef[19] * t28);
	double t184 = Bm * Bm;
	double t185 = t184 * t184;
	double t186 = 0.1e1 / t185;
	double t194 = 0.1e1 / t184 / Bm;
	double t202 = 0.1e1 / t185 / Bm;
	double t258 = 0.1e1 / t184;
	complex<double> t353 = Nm[0] * (t186 * coefm[16] * t28 + t186 * coefm[17] * t15 + t194 * coefm[1] * t7 + t186 * coefm[0] * t22 + t202 * coefm[37] * t36 + t202 * coefm[38] * t24) + Nm[5] * (t186 * coefm[34] * t22 + t186 * coefm[35] * t15 + t194 * coefm[15] * t7 + t186 * coefm[36] * t28) + Nm[8] * t194 * t7 * coefm[39] + Nm[6] * t194 * t7 * coefm[18] + Nm[1] * (t186 * coefm[32] * t15 + t186 * coefm[33] * t22 + t202 * coefm[4] * t24 + t202 * coefm[6] * t36 + t194 * coefm[19] * t7 + t186 * coefm[42] * t28) + Nm[4] * (t258 * coefm[24] * t36 + t258 * coefm[25] * t24) + Nm[11] * (t258 * coefm[20] * t22 + t258 * coefm[21] * t15 + t258 * coefm[23] * t28) + Nm[3] * (t258 * coefm[7] * t22 + t258 * coefm[13] * t15 + t194 * coefm[30] * t36 + t194 * coefm[31] * t24 + t258 * coefm[3] * t28) + Nm[2] * (t194 * coefm[8] * t15 + t194 * coefm[9] * t22 + t186 * coefm[28] * t36 + t186 * coefm[29] * t24 + t258 * coefm[41] * t7 + t194 * coefm[10] * t28) + Nm[9] * (t194 * coefm[26] * t22 + t194 * coefm[27] * t15 + t258 * coefm[5] * t7 + t194 * coefm[22] * t28) + Nm[7] * (t186 * coefm[11] * t15 + t186 * coefm[14] * t22 + t194 * coefm[40] * t7 + t186 * coefm[2] * t28) + Nm[10] * t258 * t7 * coefm[12];
	//
	X = t353 + t182;
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> X_function_f3_f3
// ***********************************************************************

complex<double> X_function_f3_f3 (double theta, double Psi, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[])
{
	complex<double> X;
	//
	double t2 = B * B;
	double t4 = 0.1e1 / t2 / B;
	double t7 = cos(Psi);
	double t9 = t2 * t2;
	double t10 = 0.1e1 / t9;
	double t13 = sin(Psi);
	double t14 = t7 * t13;
	double t17 = 0.1e1 / t9 / B;
	double t20 = t7 * t7;
	double t21 = cos(theta);
	double t22 = t20 * t21;
	double t23 = t22 * t13;
	double t30 = sin(theta);
	double t31 = t20 * t30;
	double t36 = t20 * t13 * t30;
	double t61 = 0.1e1 / t2;
	complex<double> t182 = N[0] * (t4 * coef[21] * t7 + t10 * coef[1] * t14 + t17 * coef[6] * t23 + t10 * coef[3] * t22 + t10 * coef[4] * t31 + t17 * coef[7] * t36) + N[5] * (t10 * coef[8] * t22 + t10 * coef[35] * t31 + t4 * coef[2] * t7 + t10 * coef[36] * t14) + N[6] * t4 * t7 * coef[5] + N[11] * (t61 * coef[9] * t22 + t61 * coef[10] * t31 + t61 * coef[12] * t14) + N[4] * (t61 * coef[13] * t23 + t61 * coef[14] * t36) + N[1] * (t10 * coef[37] * t31 + t10 * coef[38] * t22 + t17 * coef[27] * t23 + t17 * coef[28] * t36 + t4 * coef[0] * t7 + t10 * coef[39] * t14) + N[8] * t4 * t7 * coef[40] + N[9] * (t4 * coef[15] * t22 + t4 * coef[16] * t31 + t61 * coef[24] * t7 + t4 * coef[11] * t14) + N[3] * (t61 * coef[29] * t22 + t61 * coef[31] * t31 + t4 * coef[19] * t23 + t4 * coef[20] * t36 + t61 * coef[25] * t14) + N[7] * (t10 * coef[32] * t22 + t10 * coef[33] * t31 + t4 * coef[41] * t7 + t10 * coef[26] * t14) + N[2] * (t4 * coef[30] * t31 + t4 * coef[34] * t22 + t10 * coef[17] * t23 + t10 * coef[18] * t36 + t61 * coef[42] * t7 + t4 * coef[22] * t14) + N[10] * t61 * t7 * coef[23];
	double t184 = Bm * Bm;
	double t185 = t184 * t184;
	double t187 = 0.1e1 / t185 / Bm;
	double t194 = 0.1e1 / t185;
	double t199 = 0.1e1 / t184 / Bm;
	double t227 = 0.1e1 / t184;
	complex<double> t353 = Nm[0] * (t187 * coefm[32] * t36 + t187 * coefm[33] * t23 + t194 * coefm[0] * t22 + t199 * coefm[31] * t7 + t194 * coefm[15] * t14 + t194 * coefm[14] * t31) + Nm[5] * (t194 * coefm[36] * t31 + t194 * coefm[37] * t22 + t199 * coefm[28] * t7 + t194 * coefm[34] * t14) + Nm[10] * t227 * t7 * coefm[5] + Nm[9] * (t199 * coefm[22] * t22 + t199 * coefm[23] * t31 + t227 * coefm[9] * t7 + t199 * coefm[18] * t14) + Nm[7] * (t194 * coefm[6] * t22 + t194 * coefm[8] * t31 + t199 * coefm[40] * t7 + t194 * coefm[3] * t14) + Nm[8] * t199 * t7 * coefm[41] + Nm[6] * t199 * t7 * coefm[30] + Nm[4] * (t227 * coefm[20] * t36 + t227 * coefm[21] * t23) + Nm[1] * (t194 * coefm[35] * t31 + t194 * coefm[42] * t22 + t187 * coefm[2] * t23 + t187 * coefm[4] * t36 + t199 * coefm[29] * t7 + t194 * coefm[38] * t14) + Nm[11] * (t227 * coefm[16] * t22 + t227 * coefm[17] * t31 + t227 * coefm[19] * t14) + Nm[3] * (t227 * coefm[7] * t22 + t227 * coefm[12] * t31 + t199 * coefm[26] * t36 + t199 * coefm[27] * t23 + t227 * coefm[1] * t14) + Nm[2] * (t199 * coefm[11] * t22 + t199 * coefm[13] * t31 + t194 * coefm[24] * t36 + t194 * coefm[25] * t23 + t227 * coefm[39] * t7 + t199 * coefm[10] * t14);
	//
	X = t353 + t182;
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF void X_function_pre
// ***********************************************************************

void X_function_pre (double theta, double Psi, double tPsiA, double tPsiB, double PsiA, double PsiB, double B, double Bm, complex<double> N[], complex<double> Nm[], const double ko)
{
	double  D, D1, D2;
	//
	complex<double> a  = Iunit * ko * B;
	complex<double> am = Iunit * ko * Bm;
	//
	if (Psi >= PsiB && Psi >= PsiA)
	{
		D  = sqrt(3.0) / sin(Psi);
		X2(D, a, N);
		X2 (D, am, Nm);
	}
	else if (theta >= M_PI_2)
	{
		D  = sqrt(3.0) / (cos(Psi) * tPsiA);
		X2 (D, a, N);
		X2 (D, am, Nm);
	}
	else if (Psi >= PsiA)
	{
		D1 = 2.0 * sqrt(3.0) / (cos(Psi) * ( tPsiB + tan(Psi) )  );
		D2 = sin(Psi) / sqrt(3.0);
		X1(Psi, D1, D2, tPsiB, a, N);
		X1(Psi, D1, D2, tPsiB, am, Nm);
	}
	else
	{
		D1 = sqrt(3.0) / (cos(Psi) * sin(theta) );
		D2 = ( cos(Psi) * tPsiA ) / sqrt(3.0);
		X1(Psi, D1, D2, tPsiB, a, N);
		X1(Psi, D1, D2, tPsiB, am, Nm);
	}

}




