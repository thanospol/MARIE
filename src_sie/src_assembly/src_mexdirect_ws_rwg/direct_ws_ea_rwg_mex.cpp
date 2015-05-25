#include "mex.h"
#include "direct_ws_ea_rwg.h"


/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks you out of
     the MEX-file) */
  if(nrhs!=11) 
    mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumInputs",
            "11 inputs required.");
  if(nlhs!=1) 
    mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumOutputs",
            "1 output required.");   
  
  // INPUT

  //  create a pointer to the input vectors 
  double* r1 = mxGetPr(prhs[0]);
  double* r2 = mxGetPr(prhs[1]);
  double* r3 = mxGetPr(prhs[2]);
  double* r4 = mxGetPr(prhs[3]);
      
  //  get the scalar input ko 
  double ko = mxGetScalar(prhs[4]);
  
  //  get the scalar input Np_1D 
  int N_theta = mxGetScalar(prhs[5]);
  int N_psi   = mxGetScalar(prhs[6]);
    
  //  create a pointer to w,z 
  double* z_theta = mxGetPr(prhs[7]);
  double* w_theta = mxGetPr(prhs[8]);
  double* z_psi   = mxGetPr(prhs[9]);
  double* w_psi   = mxGetPr(prhs[10]);
  
//   mexPrintf("z_psi =  (%f,%f,%f) \n", z_psi[7],z_psi[8],z_psi[9]);
  
  // OUTPUT
  //  set the output pointer to the output matrix 
  plhs[0] = mxCreateDoubleMatrix(9, 1, mxCOMPLEX);    
  
  //  create a C pointer to a copy of the output matrix 
  double* I_DEr = mxGetPr(plhs[0]);
  double* I_DEi = mxGetPi(plhs[0]);
  
  /*  call the C subroutine */
  complex<double> I_DE[9];
  
   direct_ws_ea_rwg( r1, r2, r3, r4, ko, N_theta, N_psi, z_theta, w_theta, z_psi, w_psi, I_DE );
 for(int i=0; i < 9; i++){
  I_DEr[i] = real(I_DE[i]);
  I_DEi[i] = imag(I_DE[i]);
  }
}
