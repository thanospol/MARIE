// #include <iostream>
// #include <iomanip>
// #include <fstream>
#include <complex>
#include <math.h>
// #include <time.h>
// #include "globals.h"

using namespace std;


void Coupling ( const double r1[3], const double r2[3], const double r3[3], const double r_obs[3], const double ko, const int Np_2D, const double Z1[], const double Z2[], const double Z3[], const double wp[], complex<double> GJ_1[3], complex<double> GJ_2[3], complex<double> GJ_3[3]  )
{
    complex<double> Iunit = complex<double>( 0.0 , 1.0 );
    //double M_PI = 3.14159265358979323846264338328 ;
    //double eo	     =   8.85400e-12;				// free space electric permitivity
    //double mo	     =   4.0 * M_PI * 1.0e-7;		// free space magnetic permeability
    //complex<double> omega = ko / sqrt(eo*mo);

    
    complex<double> mul_ct ;
    
    double r_src[3], L1[3], L2[3], L3[3], rho1[3], rho2[3], rho3[3];
    double L1_length, L2_length, L3_length;
    double R_vec[3];
    double R;
    complex<double> kappa, P, Q;
    complex<double> Gxx, Gxy, Gxz, Gyy, Gyz, Gzz;
    complex<double> Kernel_1[3], Kernel_2[3], Kernel_3[3];
    //
    complex<double> Integral_1[3], Integral_2[3], Integral_3[3];
    //
    for (int i = 0; i < 3; i++)
    {
        L1[i] = r2[i] - r3[i];
        L2[i] = r3[i] - r1[i];
        L3[i] = r1[i] - r2[i];
        //
        Integral_1[i] = 0.0;
        Integral_2[i] = 0.0;
        Integral_3[i] = 0.0;
    }
    //
    L1_length = sqrt(L1[0]*L1[0] + L1[1]*L1[1] + L1[2]*L1[2]);
    L2_length = sqrt(L2[0]*L2[0] + L2[1]*L2[1] + L2[2]*L2[2]);
    L3_length = sqrt(L3[0]*L3[0] + L3[1]*L3[1] + L3[2]*L3[2]);
    //
    for (int i = 0; i < Np_2D; i++)
    {
        for (int iv = 0; iv < 3; iv++)
        {
            r_src[iv] = r1[iv] * Z1[i] + r2[iv] * Z2[i] + r3[iv] * Z3[i];
            R_vec[iv] = r_obs[iv] - r_src[iv];
            //
            rho1[iv] = r_src[iv] - r1[iv];
            rho2[iv] = r_src[iv] - r2[iv];
            rho3[iv] = r_src[iv] - r3[iv];
        }
        R = sqrt(R_vec[0]*R_vec[0] + R_vec[1]*R_vec[1] + R_vec[2]*R_vec[2]);
        //
        kappa =  exp(-Iunit * ko * R) /  R / R / R;
        P = (ko*ko*R*R - 3.0*Iunit*ko*R - 3.0) / R / R;
        Q = (ko*ko*R*R - Iunit*ko*R - 1.0) / P;
        
        mul_ct = kappa * P;
        //
        Gxx = R_vec[0] * R_vec[0] - Q;
        Gxy = R_vec[0] * R_vec[1];
        Gxz = R_vec[0] * R_vec[2];
        //
        Gyy = R_vec[1] * R_vec[1] - Q;
        Gyz = R_vec[1] * R_vec[2];
        //
        Gzz = R_vec[2] * R_vec[2] - Q;
        //
        Kernel_1[0] = Gxx * rho1[0] +  Gxy * rho1[1] + Gxz * rho1[2];
        Kernel_1[1] = Gxy * rho1[0] +  Gyy * rho1[1] + Gyz * rho1[2];
        Kernel_1[2] = Gxz * rho1[0] +  Gyz * rho1[1] + Gzz * rho1[2];

        Kernel_2[0] = Gxx * rho2[0] +  Gxy * rho2[1] + Gxz * rho2[2];
        Kernel_2[1] = Gxy * rho2[0] +  Gyy * rho2[1] + Gyz * rho2[2];
        Kernel_2[2] = Gxz * rho2[0] +  Gyz * rho2[1] + Gzz * rho2[2];

        Kernel_3[0] = Gxx * rho3[0] +  Gxy * rho3[1] + Gxz * rho3[2];
        Kernel_3[1] = Gxy * rho3[0] +  Gyy * rho3[1] + Gyz * rho3[2];
        Kernel_3[2] = Gxz * rho3[0] +  Gyz * rho3[1] + Gzz * rho3[2];
        //
        Integral_1[0] += wp[i] * mul_ct * Kernel_1[0];
        Integral_1[1] += wp[i] * mul_ct * Kernel_1[1];
        Integral_1[2] += wp[i] * mul_ct * Kernel_1[2];

        Integral_2[0] += wp[i] * mul_ct * Kernel_2[0];
        Integral_2[1] += wp[i] * mul_ct * Kernel_2[1];
        Integral_2[2] += wp[i] * mul_ct * Kernel_2[2];

        Integral_3[0] += wp[i] * mul_ct * Kernel_3[0];
        Integral_3[1] += wp[i] * mul_ct * Kernel_3[1];
        Integral_3[2] += wp[i] * mul_ct * Kernel_3[2];
    }

    
    //
    for (int i = 0; i < 3; i++)
    {
        GJ_1[i] = L1_length  * Integral_1[i];
        GJ_2[i] = L2_length  * Integral_2[i];
        GJ_3[i] = L3_length  * Integral_3[i];
    }

}
