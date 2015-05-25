/**************************************************************************************************************************
           
		           globas.h

Global constants, definitions & function^s declaration.

  Licensing: This code is distributed under the GNU LGPL license. 

  Modified:  04 February 2011

  Author:    Athanasios Polimeridis

  Reference

  A. G. Polimeridis, and T. V. Yioultsis, “On the direct evaluation of weakly singular
  integrals in Galerkin mixed potential integral equation formulations,” IEEE Trans.
  Antennas Propag., vol. 56, no. 9, pp. 3011-3019, Sep. 2008.

  A. G. Polimeridis, and J. R. Mosig, “Complete semi-analytical treatment of weakly
  singular integrals on planar triangles via the direct evaluation method,” Int. J.
  Numerical Methods Eng., vol. 83, pp. 1625-1650, 2010.

**************************************************************************************************************************/

#if !defined _GLOBALS_H_
#define _GLOBALS_H_

#include <math.h>
#include <complex>

using namespace std;

// typedef long double MYFLOAT;

// *********************************************
//			DECLARATION OF FUNCTIONS
// *********************************************

complex<double> directfn_ws_st_const (const double r1[], const double r2[], const double r3[], const double ko, const int Np_1D );

void PSI_limits ( int argument, double *psi_A, double *psi_B );
void ETA_limits ( int argument, double PSI, double *ETA_A, double *ETA_B );

complex<double> N_functions ( const double r1[], const double r2[], const double r3[], const double ko, int basis_testing, int argument, double PSI, double ETA, int kk, const int Np_1D, const double z_LAMBDA[], const double w_LAMBDA[] );
complex<double> A_functions ( const double r1[], const double r2[], const double r3[], const double ko, int basis_testing, double PSI, double ETA, double LAMBDA, int argument, const int Np_1D, const double z_RHO[], const double w_RHO[] );

void Subtriangles(double ETA_sub, double XI_sub, int argument, double *ETA, double *XI);

void Simplex(double *xi_1p, double *xi_2p, double *xi_3p, double *xi_1q, double *xi_2q, double *xi_3q, double ETA_p, double XI_p, double ETA_q, double XI_q);

complex<double> Kernel (const double r1[], const double r2[], const double r3[], const double ko, int argument, double xi_1p, double xi_2p, double xi_3p, double xi_1q, double xi_2q, double xi_3q );


void GL_1D ( int n, double x[], double w[] );

//  INPUT

/*const double r1[] = {0.0 , 0.0 , 0.5};
const double r2[] = {0.0 , 0.0 , 0.0};
const double r3[] = {0.5 , 0.0 , 0.0};

const double ko = double(1.0);*/

// ********
//  Macros
// ********

#define round(x)	        ( (int)((x)+ 0.5) )
#define RoundToZero(x)		( fabs(x) > 1.0e-12 ? (x) : 0.0 )

// **************************************
//			Inline functions
// **************************************

inline
        double vector_dot(double x[], double y[]) {
        return x[0]*y[0]+x[1]*y[1]+x[2]*y[2];
}
//
inline
        void vector_cross(double x[], double y[], double z[])
    {
        z[0] = x[1] * y[2] - x[2] * y[1];
        z[1] = x[2] * y[0] - x[0] * y[2];
        z[2] = x[0] * y[1] - x[1] * y[0];
    }

// **************************************
//			Mathematical Constants
// **************************************

#ifndef M_PI
#define M_PI       3.14159265358979323846264338328      /* pi */
#endif

#ifndef M_PI_2
#define M_PI_2     1.57079632679489661923132169164      /* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4     0.78539816339744830961566084582     /* pi/4 */
#endif

#ifndef M_SQRTPI
#define M_SQRTPI   1.77245385090551602729816748334      /* sqrt(pi) */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257389615890312      /* 2/sqrt(pi) */
#endif

#ifndef M_1_PI
#define M_1_PI     0.31830988618379067153776752675      /* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI     0.63661977236758134307553505349      /* 2/pi */
#endif

// **************************************
//			Physical Constants
// **************************************

const double eo	     =   8.85400e-12;			// free space electric permitivity
const double mo	     =   4.0 * M_PI * 1.0e-7;		// free space magnetic permeability
const double co             =   299792458;			// free spave light velocity
const double Zo             =   376.734;			// free space inpedance



const std::complex<double> Iunit = std::complex<double>( 0.0 , 1.0 );

#endif /* _GLOBALS_H_ */
