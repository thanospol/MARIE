#include "mex.h"
#include <math.h>
#include "matrix.h"
#include <complex>

using namespace std;

void cubature_nop ( const int Np, const double w[], const double u[], const double v[], const double r_m[3], const double r_n[3], const double R_faces[3][4][6], const double ko, const double dx, complex<double> I_mn[6][6]  );

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks you out of
     the MEX-file) */
  if(nrhs!=9) 
    mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumInputs",
            "11 inputs required.");
  if(nlhs!=1) 
    mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumOutputs",
            "1 output required.");   
  
  // INPUT
  
  //  get the scalar input Np 
  int Np = (int) mxGetScalar(prhs[0]);
  //  create a pointer to w,u,v 
  double* w = mxGetPr(prhs[1]);
  double* u = mxGetPr(prhs[2]);
  double* v = mxGetPr(prhs[3]);

  //  create a pointer to the input vectors r_m, r_n 
  double* r_m = mxGetPr(prhs[4]);
  double* r_n = mxGetPr(prhs[5]);

  // const double R_faces[3][4][6]
  int dimensions[3];
  dimensions[0] = 3;
  dimensions[1] = 4;
  dimensions[2] = 6;

  double* R_faces = mxGetPr(prhs[6]);
  
  double RR_faces[3][4][6];
  
  for(int k=0; k < 6; k++)
   {
       for(int j=0; j < 4; j++)
       {    
           for(int i=0; i < 3; i++)
           {
             RR_faces[i][j][k]  = R_faces[i + 3*j + 3*4*k] ;
           }
          
       }
   }
  
  //  get the scalar input ko 
  double ko = mxGetScalar(prhs[7]);
  
  //  get the scalar input dx 
  double dx = mxGetScalar(prhs[8]);
  
  // OUTPUT
  //  set the output pointer to the output matrix 
  plhs[0] = mxCreateDoubleMatrix(6, 6, mxCOMPLEX);    
  
  //  create a C pointer to a copy of the output matrix 
  double* I_mnr = mxGetPr(plhs[0]);
  double* I_mni = mxGetPi(plhs[0]);
  
  /*  call the C subroutine */
  complex<double> I_mn[6][6];
  
   cubature_nop ( Np, w, u, v, r_m, r_n, RR_faces, ko, dx, I_mn  );
   //
   for(int i=0; i < 6; i++)
   {
       for(int j=0; j < 6; j++)
       {        
          I_mnr[i + 6*j] = real(I_mn[i][j]);
          I_mni[i + 6*j] = imag(I_mn[i][j]);
       }
   }
   
}
