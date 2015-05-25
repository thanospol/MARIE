#include "mex.h"
#include "direct_ws_va_const.h"

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks you out of
     the MEX-file) */
  if(nrhs!=15) 
    mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumInputs",
            "15 inputs required.");
  if(nlhs!=1) 
    mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumOutputs",
            "1 output required.");   
  
  // INPUT

  //  create a pointer to the input vectors 
  double* r1 = mxGetPr(prhs[0]);
  double* r2 = mxGetPr(prhs[1]);
  double* r3 = mxGetPr(prhs[2]);
  double* r4 = mxGetPr(prhs[3]);
  double* r5 = mxGetPr(prhs[4]);
    
  
  //  get the scalar input ko 
  double ko = mxGetScalar(prhs[5]);
  
  //  get the scalar input Np_1D 
  int N_theta_p = (int) mxGetScalar(prhs[6]);
  int N_theta_q = (int) mxGetScalar(prhs[7]);
  int N_psi     = (int) mxGetScalar(prhs[8]);
    
  //  create a pointer to w,z 
  double* z_theta_p = mxGetPr(prhs[9]);
  double* w_theta_p = mxGetPr(prhs[10]);
  double* z_theta_q = mxGetPr(prhs[11]);
  double* w_theta_q = mxGetPr(prhs[12]);
  double* z_psi     = mxGetPr(prhs[13]);
  double* w_psi     = mxGetPr(prhs[14]);
  
//   mexPrintf("z_psi =  (%f,%f,%f) \n", z_psi[7],z_psi[8],z_psi[9]);
  
  // OUTPUT
  //  set the output pointer to the output matrix 
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxCOMPLEX);  
  
  //  create a C pointer to a copy of the output matrix 
  double* I_DEr = mxGetPr(plhs[0]);
  double* I_DEi = mxGetPi(plhs[0]);
  
  /*  call the C subroutine */
  complex<double> I_DE = direct_ws_va_const( r1, r2, r3, r4, r5, ko, N_theta_p, N_theta_q, N_psi, z_theta_p, w_theta_p, z_theta_q, w_theta_q, z_psi, w_psi );
  I_DEr[0] = real(I_DE);
  I_DEi[0] = imag(I_DE);
}
