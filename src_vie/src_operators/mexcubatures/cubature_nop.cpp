#include <complex>
#include <math.h>
//#include "globals.h"

using namespace std;


void cubature_nop ( const int Np, const double w[], const double u[], const double v[], const double r_m[3], const double r_n[3], const double R_faces[3][4][6], const double ko, const double dx, complex<double> I_mn[6][6]  )
{
	complex<double> Iunit = complex<double>( 0.0 , 1.0 );
	//
    int i_obs, i_src;
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
	//
    double Wt;
    complex<double> Integral;
 // int index_obs[21] = {1,1,1,1,1,1,2,2,2,2,2,3,3,3,3,4,4,4,5,5,6};
	int index_obs[21] = {0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,3,3,3,4,4,5};

  //int index_src[21] = {1,2,3,4,5,6,1,3,4,5,6,3,4,5,6,3,5,6,5,6,5};
	int index_src[21] = {0,1,2,3,4,5,0,2,3,4,5,2,3,4,5,2,4,5,4,5,4};
	// ********************************************
    //		 Collect all 21 unique interactions
    // ********************************************
    for (int i = 0; i < 21; i++)
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
        Integral = 0.0;
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
                    Rmn[j] = rrq[j] - rrp[j];

                }
                //
                //Rpq = sqrt(vector_dot(Rmn, Rmn));
				Rpq = sqrt(Rmn[0]*Rmn[0]+Rmn[1]*Rmn[1]+Rmn[2]*Rmn[2]);
                //
                G = exp(-Iunit * ko * Rpq) / Rpq;
                Integral += Wt * G;

            }

        }
		//
        I_mn[i_obs][i_src] = Integral;	
    }
    // The non-unique elements
    I_mn[1][1] = I_mn[0][0];
    //
    I_mn[2][0] = I_mn[1][3];
    I_mn[2][1] = I_mn[0][3];
    //
    I_mn[3][0] = I_mn[1][2];
    I_mn[3][1] = I_mn[0][2];
    I_mn[3][3] = I_mn[2][2];
    //
    I_mn[4][0] = I_mn[1][5];
    I_mn[4][1] = I_mn[0][5];
    I_mn[4][2] = I_mn[3][5];
    I_mn[4][3] = I_mn[2][5];
    //
    I_mn[5][0] = I_mn[1][4];
    I_mn[5][1] = I_mn[0][4];
    I_mn[5][2] = I_mn[3][4];
    I_mn[5][3] = I_mn[2][4];
    I_mn[5][5] = I_mn[4][4];
	//
	double pow4 = dx*dx*dx*dx;

	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			I_mn[i][j] = pow4 * I_mn[i][j];
		}
	}
}
