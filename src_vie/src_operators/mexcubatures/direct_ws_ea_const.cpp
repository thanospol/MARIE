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

#include "direct_ws_ea_const.h"

// ***********************************************************************
//			IMPLEMENTATION OF void DIRECT
// ***********************************************************************

complex<double> direct_ws_ea_const (const double r1[],const double r2[],const double r3[],const double r4[], const double ko, const int N_theta, const int N_psi, const double w_theta[], const double z_theta[], const double w_psi[], const double z_psi[] )
{
	// ************************************************
	//			DECLARATION OF KEY VARIABLES
	// ************************************************

	// 1. Various

	double alpha[3], beta[3], gamma[3], r21[3], r31[3], r32[3], r41[3], r42[3];
	//
	complex<double> Ising_const;
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
	complex<double> I_theta_const;
	complex<double> I_psi_const;
	//
	complex<double> X_const;
	//
	complex<double> I_const[6];
	//
	complex<double> Iconst;

	// 2. Coefficients' parameters

	complex<double> coef_const[3],  coefm_const[3];

	// 3.

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
    // Compute area
    double Ap_[3];
    vector_cross(r21, r42, Ap_);
    double Ap = 1.0 / 2.0 * sqrt(vector_dot(Ap_, Ap_));
    double Jp = Ap / sqrt(3.0);
    //
    double Aq_[3];
    vector_cross(r41, r42, Aq_);
    double Aq = 1.0 / 2.0 * sqrt(vector_dot(Aq_, Aq_));
    double Jq = Aq / sqrt(3.0);
    double J = Jp * Jq;
	// Get the coefficients
	coefficients_const ( r1,  r2,  r3,  r4,  ko,  coef_const,  coefm_const );
     // Initialization of I_
	 for ( int im = 0; im <  6; im++ )
	 {
		 I_const[im] = 0.0;
	 }
	 //
	 for ( int m = 1; m <  7; m++ )
	 {
		 I_theta_const = 0.0;
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
				 // Integrate
				 WPSI = w_psi[n_psi];
				 //
				 I_psi_const +=  WPSI * X_const;
			 }//end for ( int n_psi = 0 ; n_psi <  N_psi ; n_psi++ )
			 //
			 I_psi_const *=  J_psi;
			 //
			 WTHETA = w_theta[n_theta];
			 //
			 I_theta_const +=  WTHETA * I_psi_const;
		 } //end for ( int n_theta = 0 ; n_theta <  N_theta ; n_theta++ )
		 //
		 I_const[m-1] = J_theta * I_theta_const;
	 } //end for ( int m = 1; m <  7; m++ )
	 //
	 Iconst = I_const[0] + I_const[1] + I_const[2] + I_const[3] + I_const[4] + I_const[5];
	 //
	 // FINAL OUTPUT
	 complex<double> I_DE = J * Iconst;
    return I_DE;
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