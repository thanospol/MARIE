#include "mex.h"
#include "direct_ws_st_rwg.h"


/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks you out of
     the MEX-file) */
  if(nrhs!=7) 
    mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumInputs",
            "7 inputs required.");
  if(nlhs!=1) 
    mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumOutputs",
            "1 output required.");   
  
  // INPUT

  //  create a pointer to the input vectors 
  double* r1 = mxGetPr(prhs[0]);
  double* r2 = mxGetPr(prhs[1]);
  double* r3 = mxGetPr(prhs[2]);
    
  //  get the scalar input ko 
  double ko = mxGetScalar(prhs[3]);
  
  //  get the scalar input Np_1D 
  int Np_1D = mxGetScalar(prhs[4]);
    
  //  create a pointer to w,z 
  double* z = mxGetPr(prhs[5]);
  double* w = mxGetPr(prhs[6]);
  
//   mexPrintf("z =  (%f,%f,%f) \n", z[0],z[1],z[2]);
  
  // OUTPUT
  //  set the output pointer to the output matrix 
  plhs[0] = mxCreateDoubleMatrix(9, 1, mxCOMPLEX);  
  
  //  create a C pointer to a copy of the output matrix 
  double* I_DEr = mxGetPr(plhs[0]);
  double* I_DEi = mxGetPi(plhs[0]);
  
  /*  call the C subroutine */
  complex<double> I_DE[9];
  
  direct_ws_st_rwg( r1, r2, r3, ko, Np_1D, z, w, I_DE );
 for(int i=0; i < 9; i++){
  I_DEr[i] = real(I_DE[i]);
  I_DEi[i] = imag(I_DE[i]);
  }
}
