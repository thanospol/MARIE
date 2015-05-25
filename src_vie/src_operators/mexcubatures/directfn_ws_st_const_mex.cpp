#include "mex.h"
#include "directfn_ws_st_const.h"


/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks you out of
     the MEX-file) */
  if(nrhs!=5) 
    mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumInputs",
            "5 inputs required.");
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
  int Np_1D = (int) mxGetScalar(prhs[4]);
    
//   //  create a pointer to w,z 
//   double* z = mxGetPr(prhs[5]);
//   double* w = mxGetPr(prhs[6]);
  
//   mexPrintf("z =  (%f,%f,%f) \n", z[0],z[1],z[2]);
  
  // OUTPUT
  //  set the output pointer to the output matrix 
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxCOMPLEX);  
  
  //  create a C pointer to a copy of the output matrix 
  double* I_DEr = mxGetPr(plhs[0]);
  double* I_DEi = mxGetPi(plhs[0]);
  
  /*  call the C subroutine */
  complex<double> I_DE = directfn_ws_st_const ( r1, r2, r3, ko, Np_1D );
  I_DEr[0] = real(I_DE);
  I_DEi[0] = imag(I_DE);
}
