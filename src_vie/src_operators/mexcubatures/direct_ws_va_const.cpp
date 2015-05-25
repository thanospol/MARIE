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

#include "direct_ws_va_const.h"
using namespace std;

// ***********************************************************************
//			IMPLEMENTATION OF void DIRECT
// ***********************************************************************

complex<double> direct_ws_va_const (const double r1[],const double r2[],const double r3[],const double r4[], const double r5[], const double ko, const int N_theta_p, const int N_theta_q, const int N_psi, const double w_theta_p[], const double z_theta_p[], const double w_theta_q[], const double z_theta_q[], const double w_psi[], const double z_psi[] )
{
	// ************************************************
	//			DECLARATION OF KEY VARIABLES
	// ************************************************

	// 1. Various

	double r21[3], r31[3], r32[3], r41[3], r54[3], r51[3];
	//
	complex<double> Ising_const;
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

    complex<double> a1, a2;
	complex<double> K_1, K_2;

	complex<double> Omega_const_1;
	complex<double> Omega_const_2;
	//
	double WPSI, WTHETA_p, WTHETA_q;
	//
	complex<double> I_theta_p_const; 
	complex<double> I_theta_q_const; 
	complex<double> I_psi_const; 
	//
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
    // Compute area
    double Ap_[3];
    vector_cross(r21, r31, Ap_);
    double Ap = 1.0 / 2.0 * sqrt(vector_dot(Ap_, Ap_));
    double Jp = Ap / sqrt(3.0);
    //
    double Aq_[3];
    vector_cross(r41, r51, Aq_);
    double Aq = 1.0 / 2.0 * sqrt(vector_dot(Aq_, Aq_));
    double Jq = Aq / sqrt(3.0);
    double J = Jp * Jq;
	// Get the coefficients
// 	coefficients_const ( r1,  r2,  r3,  r4, r5,  ko,  coef_const );
    double coef_const = 1.0;
     // Initialization of I_theta_p
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
                // K_function defined as inline in .h file
                a1 = Iunit * ko * B_1;
                K_1 = K_function(a1, L1);
                                
                a2 = Iunit * ko * B_2;
                K_2 = K_function(a2, L2);
                //
                Omega_const_1 = K_1 * (sin(PSI_1) * cos(PSI_1) * coef_const) / B_1;
                Omega_const_2 = K_2 * (sin(PSI_2) * cos(PSI_2) * coef_const) / B_2;
                //
                WPSI = w_psi[n_psi];
                //
                I_psi_const +=  WPSI * (J_psi_1 * Omega_const_1 + J_psi_2 * Omega_const_2);
            }//end for ( int n_psi = 0 ; n_psi <  N_psi ; n_psi++ )
            WTHETA_q = w_theta_q[n_theta_q];
            //
            I_theta_q_const   +=  WTHETA_q * I_psi_const;
            
        } //end for ( int n_theta_q = 0 ; n_theta_q <  N_theta_q ; n_theta_q++ )
        //
        I_theta_q_const   *=  J_theta_q ;
        //
        WTHETA_p = w_theta_p[n_theta_p];
        //
        I_theta_p_const  += WTHETA_p * I_theta_q_const;
        
     } //end for ( int n_theta_p = 0 ; n_theta_p <  N_theta_p ; n_theta_p++ )
        //
        I_theta_p_const   *=  J_theta_p ;
        // FINAL OUTPUT
        complex<double> I_DE = J * I_theta_p_const;
        return I_DE;
}