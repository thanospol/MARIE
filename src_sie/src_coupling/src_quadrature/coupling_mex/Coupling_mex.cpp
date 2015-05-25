#include "mex.h"
#include <math.h>
#include "matrix.h"
#include <complex>

using namespace std;

void Coupling ( const double r1[3], const double r2[3], const double r3[3], const double r_obs[3], const double ko, const int Np_2D, const double Z1[], const double Z2[], const double Z3[], const double wp[], complex<double> GJ_1[3], complex<double> GJ_2[3], complex<double> GJ_3[3]  );

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks you out of
     the MEX-file) */
  if(nrhs!=10) 
    mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumInputs",
            "10 inputs required.");
  if(nlhs!=3) 
    mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumOutputs",
            "3 output required.");   
  
  // INPUT
  
  //  create a pointer to the input vectors r_m, r_n 
  double* r1    = mxGetPr(prhs[0]);
  double* r2    = mxGetPr(prhs[1]);
  double* r3    = mxGetPr(prhs[2]);
  double* r_obs = mxGetPr(prhs[3]);
  
  //  get the scalar input ko 
  double ko = mxGetScalar(prhs[4]);
  
  //  get the scalar input Np 
  int Np_2D = (int) mxGetScalar(prhs[5]);
  
  //  create a pointer to w,u,v 
  double* Z1 = mxGetPr(prhs[6]);
  double* Z2 = mxGetPr(prhs[7]);
  double* Z3 = mxGetPr(prhs[8]);
  double* wp = mxGetPr(prhs[9]);
  
  
  // OUTPUT
  //  set the output pointer to the output matrix 
  plhs[0] = mxCreateDoubleMatrix(3, 1, mxCOMPLEX);    
  plhs[1] = mxCreateDoubleMatrix(3, 1, mxCOMPLEX);
  plhs[2] = mxCreateDoubleMatrix(3, 1, mxCOMPLEX);
  
  //  create a C pointer to a copy of the output matrix 
  double* GJ_1r = mxGetPr(plhs[0]);
  double* GJ_1i = mxGetPi(plhs[0]);
  
  double* GJ_2r = mxGetPr(plhs[1]);
  double* GJ_2i = mxGetPi(plhs[1]);
  
  double* GJ_3r = mxGetPr(plhs[2]);
  double* GJ_3i = mxGetPi(plhs[2]);
  
  /*  call the C subroutine */
  complex<double> GJ_1[3], GJ_2[3], GJ_3[3];
  
   Coupling ( r1, r2, r3, r_obs, ko, Np_2D, Z1, Z2, Z3, wp, GJ_1, GJ_2, GJ_3  );
   //
   for(int i=0; i < 3; i++)
   {  
      GJ_1r[i] = real(GJ_1[i]);
      GJ_1i[i] = imag(GJ_1[i]);
      
      GJ_2r[i] = real(GJ_2[i]);
      GJ_2i[i] = imag(GJ_2[i]);
      
      GJ_3r[i] = real(GJ_3[i]);
      GJ_3i[i] = imag(GJ_3[i]);
   }
   
}
