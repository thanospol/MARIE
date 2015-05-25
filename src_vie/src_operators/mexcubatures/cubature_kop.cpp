#include <complex>
#include <math.h>
//#include "globals.h"

using namespace std;


void cubature_kop ( const int Np, const double w[], const double u[], const double v[], const double r_m[3], const double r_n[3], const double R_faces[3][4][6], const double ko, complex<double> IK_mn[3]  )
{
	complex<double> Iunit = complex<double>( 0.0 , 1.0 );
	//
    int i_obs, i_src;
    double sign_ind;
    double rp1[3], rp2[3], rp4[3];
    double a_p[3], b_p[3];
    double rq1[3], rq2[3], rq4[3];
    double a_q[3], b_q[3];
	//
    double rrp[3], rrq[3];
    double Rmn[3];
    double Rpq;
    
	//
    complex<double> G;
    complex<double> w_G;
	//
    double Wt;
    complex<double> Integral_x, Integral_y, Integral_z;
    // index_obs      = [1  1  2  2  3  3  4  4  5  5  6  6]';
	int index_obs[12] = {0,0,1,1,2,2,3,3,4,4,5,5};

    // index_src      = [1  2  1  2  3  4  3  4  5  6  5  6]';
	int index_src[12] = {0,1,0,1,2,3,2,3,4,5,4,5};
    
    // sign_obs_src         = [1   -1   -1   1   1   -1   -1   1   1   -1   -1   1  ]';
	double sign_obs_src[12] = {1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0};

	// ********************************************
    //		 Collect all 21 unique interactions
    // ********************************************
    complex<double> IK_mn_x = 0.0;
    complex<double> IK_mn_y = 0.0;
    complex<double> IK_mn_z = 0.0;
    for (int i = 0; i < 12; i++)
    {
        // Observation Panel
        i_obs = index_obs[i];		
        //
        for (int j = 0; j < 3; j++)
        {
            rp1[j] = R_faces[j][0][i_obs];
            rp2[j] = R_faces[j][1][i_obs];
            rp4[j] = R_faces[j][3][i_obs];
            //
            a_p[j] = rp2[j] - rp1[j];
            b_p[j] = rp4[j] - rp1[j];
        }
        // Source Panel
        i_src = index_src[i];
        //
        for (int j = 0; j < 3; j++)
        {
            rq1[j] = R_faces[j][0][i_src];
            rq2[j] = R_faces[j][1][i_src];
            rq4[j] = R_faces[j][3][i_src];
            //
            a_q[j] = rq2[j] - rq1[j];
            b_q[j] = rq4[j] - rq1[j];
        }
        // ********************************************
        //		 4D QUADRATURE
        // ********************************************
        Integral_x = 0.0;
        Integral_y = 0.0;
        Integral_z = 0.0;
        for (int ip = 0; ip < Np; ip++)
        {
            for (int j = 0; j < 3; j++)
            {
                rrp[j] = rp1[j] + a_p[j]*u[ip] + b_p[j]*v[ip] + r_m[j];

            }

            for (int iq = 0; iq < Np; iq++)
            {
                Wt = w[ip] * w[iq];
                //
                for (int j = 0; j < 3; j++)
                {
                    rrq[j] = rq1[j] + a_q[j]*u[iq] + b_q[j]*v[iq] + r_n[j];
                    Rmn[j] = rrp[j] - rrq[j];

                }
                //
                //Rpq = sqrt(vector_dot(Rmn, Rmn));
				Rpq = sqrt(Rmn[0]*Rmn[0]+Rmn[1]*Rmn[1]+Rmn[2]*Rmn[2]);
                // (1-exp(-1i*ko*Rpq) .* (1 + 1i*ko*Rpq) )./ Rpq.^3;
                G = ( 1.0 - exp(-Iunit * ko * Rpq) * (1.0 + Iunit * ko * Rpq) ) / Rpq / Rpq / Rpq ;
                //
                w_G = Wt * G;
                Integral_x += w_G * Rmn[0];
                Integral_y += w_G * Rmn[1];
                Integral_z += w_G * Rmn[2];

            }

        }
        
        sign_ind = sign_obs_src[i];
        //
        IK_mn_x += sign_ind * Integral_x;
        IK_mn_y += sign_ind * Integral_y;
        IK_mn_z += sign_ind * Integral_z;	
    }   
    //
    IK_mn[0] = IK_mn_x;	
    IK_mn[1] = IK_mn_y;
    IK_mn[2] = IK_mn_z;
    
}
