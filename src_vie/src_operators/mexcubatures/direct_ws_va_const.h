/**************************************************************************************************************************

			Global constants, definitions & function^s declaration


  Licensing: This code is distributed under the GNU LGPL license. 

  Modified:  10 December 2010

  Author:    Athanasios Polimeridis
				
**************************************************************************************************************************/
#if !defined _GLOBALS_H_
#define _GLOBALS_H_

#include <math.h>
#include <complex>

using namespace std;

// *********************************************
//			DECLARATION OF FUNCTIONS
// *********************************************

complex<double> direct_ws_va_const (const double r1[],const double r2[],const double r3[],const double r4[], const double r5[], const double ko, const int N_theta_p, const int N_theta_q, const int N_psi, const double w_theta_p[], const double z_theta_p[], const double w_theta_q[], const double z_theta_q[], const double w_psi[], const double z_psi[]);

// ********
//  Macros
// ********

#define round(x)			( (int)((x)+ 0.5) ) 
#define RoundToZero(x)		( fabs(x) > 1.0e-12 ? (x) : 0.0 )

// **************************************
//			Inline functions
// **************************************

inline
	double vector_dot(double x[], double y[]) 
{
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
//
inline
        complex<double> K_function(complex<double> a, complex<double> L)
{
    return (1.0 / pow(a,3.0) ) * (2.0 - 2.0 * exp(-a * L) - 2.0 * a * L * exp(-a * L) - pow(a,2.0) * pow(L,2.0) * exp(-a * L) );
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

const double eo	     =   8.85400e-12;				// free space electric permitivity
const double mo	     =   4.0 * M_PI * 1.0e-7;		// free space magnetic permeability
const double co      =   299792458;			        // free spave light velocity
const double Zo		 =   376.734;					// free space inpedance

//const double zero			   = 1.0e-12;			// zero level for numerical comparisons
//const double convergence_level = 1.0e-10;			// zero level for convergence

const std::complex<double> Iunit = std::complex<double>( 0.0 , 1.0 );

#endif /* _GLOBALS_H_ */