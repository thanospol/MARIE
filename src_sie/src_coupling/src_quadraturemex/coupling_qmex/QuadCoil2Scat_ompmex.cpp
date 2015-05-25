#include "mex.h"
#include <math.h>
#include "matrix.h"
#include <complex>
#include <omp.h>
#include "stdint.h"

using namespace std;

void Coupling ( const double r1[3], const double r2[3], const double r3[3], const double r_obs[3], const double ko, const int Np_2D, const double Z1[], const double Z2[], const double Z3[], const double wp[], complex<double> GJ_1[3], complex<double> GJ_2[3], complex<double> GJ_3[3]  );

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{  
	// ------------------------------------------------------------------------------
	//
	//  Receives:
	//      R1 - double 3 x NE vector with the positions of the first vertex of each element
	//      R2 - double 3 x NE vector with the positions of the second vertex of each element
	//      R3 - double 3 x NE vector with the positions of the third vertex of each element
	//		NE - number of elements
	//      RO - double 3 x NO vector with the positions of the observation points
	//		NO - number of observation points
	//		IDX - 3 x NE integers with the index of each edge of each element (incidence of I)
	//		MULT - 3 x NE doubles with the multiplier (+1 or -1) of each edge of each element (contribution of I)
	//		NC - number of Coil edges
	//		K0 - wave number
	//		NP_2D - data for quadrature rule
	//		Z1 - more data for quadrature
	//		Z2 - and more
	//		Z3 - and yet some more
	//		WP - more?
	//		I - NC complex numbers with currents at each edge
	//
	//
	//  Returns:
	//      V - if I exists: No x 3 complex with each component of the field in each of the No positions
	//      V - if I does not exists: No x 3 x NC complex with the coupling matrix
	//
	//
	//	NOTE: 	all I/O are treated as vectors, and ordering in column major format
	//			for the coordinates, index and mult, the ordering is for ii point, x = R[3*ii], y = R[3*ii+1], z = R[3*ii+2]
	//
	// 	MATLAB CALL:
	//			[V] = mexQuadCoil2Scat(R1,R2,R3,NE,RO,NO,IDX,MULT,NC,K0,NP_2D,Z1,Z2,Z3,WP,I)
	//

	// ----------------------------------------------------------------------------
	// some initial declarations
	
	int outvec = 0; // indicates if output is the coupling matrix (0), or vector with fields (1), if I is given
	
	switch (nrhs){
		case(15):{
			outvec = 0;
			break;
		}
		case(16):{
			outvec = 1;
			break;
		}
        default:{
            mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumInputs", "QuadCoil2Scat_mex - Either 15 or 16 inputs!");
            break;
        }
	}
	
	// ----------------------------------------------------------------------------
	// declarations
	
	// aux vars
	int NPROC = 1;
	int64_t ii = 0, jj = 0, shift = 0;
	
	//  point to data in vectors R 
	double* R1 = mxGetPr(prhs[0]);
	double* R2 = mxGetPr(prhs[1]);
	double* R3 = mxGetPr(prhs[2]);
	int64_t NE = (int64_t) mxGetScalar(prhs[3]);
	double* RO = mxGetPr(prhs[4]);
	int64_t NO = (int64_t) mxGetScalar(prhs[5]);
	
	// pointer to the indexes
	double* IDXdouble = mxGetPr(prhs[6]); 
	
	// pointer to the multipliers
	double* MULT = mxGetPr(prhs[7]);
	
	// get the scalar input NC
	int64_t NC = (int64_t) mxGetScalar(prhs[8]);

	//  get the scalar input K0 
	double K0 = mxGetScalar(prhs[9]);

	//  get the scalar input NP 
	int NP_2D = (int) mxGetScalar(prhs[10]);

	//  create a pointer to w,u,v 
	double* Z1 = mxGetPr(prhs[11]);
	double* Z2 = mxGetPr(prhs[12]);
	double* Z3 = mxGetPr(prhs[13]);
	double* WP = mxGetPr(prhs[14]);
	
		
  	// this is for the computation of each coupling entry
    complex<double> GJ_1[3], GJ_2[3], GJ_3[3];
	double r1[3], r2[3], r3[3], r_obs[3];
	int64_t IDXlocal[3]; // index local
	double MULTlocal[3]; // local multiplier
	
	// for the input current vector
	double *Ireal, *Iimag; // initial declaration of pointer to I

	// for the output matrix V
	double *Vreal, *Vimag; // not yet known the size
	
	// check the data for the input current vector
	if (outvec > 0){
		// we generate the fields output
		Ireal = mxGetPr(prhs[15]);
        Iimag = mxGetPi(prhs[15]);
        
        // make sure that if the imaginary part of J is zero, this does not break
        if (Iimag == NULL){
            Iimag = (double*) calloc (NC,sizeof(double)); // just zero array
        }
				
		// allocate space for the NO x 3 complex vector
		// allocate MATLAB size for the outputs, and assign pointer
        shift = 3*NO;
        plhs[0] = mxCreateDoubleMatrix(shift, 1, mxCOMPLEX); // this initializes to 0
        Vreal = mxGetPr(plhs[0]);
        Vimag = mxGetPi(plhs[0]);
				
	}
	else{
		// we generate the coupling matrix
		// allocate space for the NO x 3 x Nedge complex vector
		// allocate MATLAB size for the outputs, and assign pointer
        shift = 3*NO;
        plhs[0] = mxCreateDoubleMatrix(shift, NC, mxCOMPLEX); // this initializes to 0
        Vreal = mxGetPr(plhs[0]);
        Vimag = mxGetPi(plhs[0]);

	}


    // ----------------------------------------------------------------------------
    // We are finally good to do some computations...
    //
    
    NPROC = omp_get_num_procs(); // get number of processors in machine
	omp_set_num_threads(NPROC); // set number of threads to maximum
	
	// The loop is reversed to the logical implementation to avoid data runs when parallel
	#pragma omp parallel for \
	default(shared) private(r_obs,r1,r2,r3,GJ_1,GJ_2,GJ_3,IDXlocal,MULTlocal,shift,jj,ii)
	for(jj = 0; jj < NO; jj++){ // loop on the number of observation points
	
		// get the coordinates of the observation point
		r_obs[0] = RO[3*jj]; // x coordinate
		r_obs[1] = RO[3*jj+1]; // y coordinate
		r_obs[2] = RO[3*jj+2]; // z coordinate		

		
		for(ii = 0; ii < NE; ii++){ // loop on the number of elements
	
			// get the coordinates of the element points
			r1[0] = R1[3*ii]; // x coordinate
			r1[1] = R1[3*ii+1]; // y coordinate
			r1[2] = R1[3*ii+2]; // z coordinate
			
			r2[0] = R2[3*ii]; // x coordinate
			r2[1] = R2[3*ii+1]; // y coordinate
			r2[2] = R2[3*ii+2]; // z coordinate
			
			r3[0] = R3[3*ii]; // x coordinate
			r3[1] = R3[3*ii+1]; // y coordinate
			r3[2] = R3[3*ii+2]; // z coordinate
		
			
			// Call the C-subroutine for each case
			Coupling ( r1, r2, r3, r_obs, K0, NP_2D, Z1, Z2, Z3, WP, GJ_1, GJ_2, GJ_3  );

			
			// get the index of the different edges
			IDXlocal[0] = (int64_t) IDXdouble[3*ii]; // index of first edge
			IDXlocal[1] = (int64_t) IDXdouble[3*ii+1]; // index of second edge
			IDXlocal[2] = (int64_t) IDXdouble[3*ii+2]; // index of third edge
			
			// get the positive or negative multiplier
			MULTlocal[0] = MULT[3*ii]; // multiplier of first edge
			MULTlocal[1] = MULT[3*ii+1]; // multiplier of second edge
			MULTlocal[2] = MULT[3*ii+2]; // multiplier of third edge
			
			
			// apply the product according to the output case
			
			if(outvec > 0){
				//vector output: we have to multiply by I
				
				// in matlab:
				// V(jj,:) = V(jj,:) + I(idx1)*mult(1)*GJ1.'; % first edge
				
				// note that the output is given in NO x 3 format, with column major ordering
				// all the I/O are vectors... reshape occurs outside		
				
				
				// first edge
				if (IDXlocal[0] > 0){ 
				
					// we have non external edge
					shift = IDXlocal[0] - 1; // recall matlab indexes start in 1
					
					//first component
					Vreal[jj] += MULTlocal[0]*( Ireal[shift]*real(GJ_1[0]) - Iimag[shift]*imag(GJ_1[0]) );
					Vimag[jj] += MULTlocal[0]*( Ireal[shift]*imag(GJ_1[0]) + Iimag[shift]*real(GJ_1[0]) );
					
					//second component (NO extra shift)
					Vreal[NO+jj] += MULTlocal[0]*( Ireal[shift]*real(GJ_1[1]) - Iimag[shift]*imag(GJ_1[1]) );
					Vimag[NO+jj] += MULTlocal[0]*( Ireal[shift]*imag(GJ_1[1]) + Iimag[shift]*real(GJ_1[1]) );
					
					//third component (2*NO extra shift)
					Vreal[2*NO+jj] += MULTlocal[0]*( Ireal[shift]*real(GJ_1[2]) - Iimag[shift]*imag(GJ_1[2]) );
					Vimag[2*NO+jj] += MULTlocal[0]*( Ireal[shift]*imag(GJ_1[2]) + Iimag[shift]*real(GJ_1[2]) );
					
				}
							
				// second edge
				if (IDXlocal[1] > 0){ 
				
					// we have non external edge
					shift = IDXlocal[1] - 1; // recall matlab indexes start in 1
					
					//first component
					Vreal[jj] += MULTlocal[1]*( Ireal[shift]*real(GJ_2[0]) - Iimag[shift]*imag(GJ_2[0]) );
					Vimag[jj] += MULTlocal[1]*( Ireal[shift]*imag(GJ_2[0]) + Iimag[shift]*real(GJ_2[0]) );
					
					//second component (NO extra shift)
					Vreal[NO+jj] += MULTlocal[1]*( Ireal[shift]*real(GJ_2[1]) - Iimag[shift]*imag(GJ_2[1]) );
					Vimag[NO+jj] += MULTlocal[1]*( Ireal[shift]*imag(GJ_2[1]) + Iimag[shift]*real(GJ_2[1]) );
					
					//third component (2*NO extra shift)
					Vreal[2*NO+jj] += MULTlocal[1]*( Ireal[shift]*real(GJ_2[2]) - Iimag[shift]*imag(GJ_2[2]) );
					Vimag[2*NO+jj] += MULTlocal[1]*( Ireal[shift]*imag(GJ_2[2]) + Iimag[shift]*real(GJ_2[2]) );
					
				}
				
				// third edge
				if (IDXlocal[2] > 0){ 
				
					// we have non external edge
					shift = IDXlocal[2] - 1; // recall matlab indexes start in 1
					
					//first component
					Vreal[jj] += MULTlocal[2]*( Ireal[shift]*real(GJ_3[0]) - Iimag[shift]*imag(GJ_3[0]) );
					Vimag[jj] += MULTlocal[2]*( Ireal[shift]*imag(GJ_3[0]) + Iimag[shift]*real(GJ_3[0]) );
					
					//second component (NO extra shift)
					Vreal[NO+jj] += MULTlocal[2]*( Ireal[shift]*real(GJ_3[1]) - Iimag[shift]*imag(GJ_3[1]) );
					Vimag[NO+jj] += MULTlocal[2]*( Ireal[shift]*imag(GJ_3[1]) + Iimag[shift]*real(GJ_3[1]) );
					
					//third component (2*NO extra shift)
					Vreal[2*NO+jj] += MULTlocal[2]*( Ireal[shift]*real(GJ_3[2]) - Iimag[shift]*imag(GJ_3[2]) );
					Vimag[2*NO+jj] += MULTlocal[2]*( Ireal[shift]*imag(GJ_3[2]) + Iimag[shift]*real(GJ_3[2]) );
					
				}

				
			} // end if outvec > 0
			else{
				// we compute the coupling matrix
				
				// in matlab
				// V(jj,:,idx1) = V(jj,:,idx1) + mult(1)*GJ1.'; % first edge
				
				// note that the output is given in NO x 3 x NC format, with column major ordering
				// all the I/O are vectors... reshape occurs outside
				
				// first edge
				if (IDXlocal[0] > 0){ 
				
					// we have non external edge
					
					// shift to the initial position of the edge related data (IDXLocal[0] in 3rd dimension)
					shift = (IDXlocal[0]-1)*3*NO; // shift 3*NO for each value of edge, recall: matlab indexes start at 1
					
					//first component
					Vreal[shift+jj] += MULTlocal[0]*real(GJ_1[0]);
					Vimag[shift+jj] += MULTlocal[0]*imag(GJ_1[0]);
					
					//second component (NO extra shift)
					Vreal[shift+NO+jj] += MULTlocal[0]*real(GJ_1[1]);
					Vimag[shift+NO+jj] += MULTlocal[0]*imag(GJ_1[1]);
					
					//third component (2*NO extra shift)
					Vreal[shift+2*NO+jj] += MULTlocal[0]*real(GJ_1[2]);
					Vimag[shift+2*NO+jj] += MULTlocal[0]*imag(GJ_1[2]);
					
				}
				
				// second edge
				if (IDXlocal[1] > 0){ 
				
					// we have non external edge
					
					// shift to the initial position of the edge related data (IDXLocal[1] in 3rd dimension)
					shift = (IDXlocal[1]-1)*3*NO; // shift 3*NO for each value of edge, recall: matlab indexes start at 1
					
					//first component
					Vreal[shift+jj] += MULTlocal[1]*real(GJ_2[0]);
					Vimag[shift+jj] += MULTlocal[1]*imag(GJ_2[0]);
					
					//second component (NO extra shift)
					Vreal[shift+NO+jj] += MULTlocal[1]*real(GJ_2[1]);
					Vimag[shift+NO+jj] += MULTlocal[1]*imag(GJ_2[1]);
					
					//third component (2*NO extra shift)
					Vreal[shift+2*NO+jj] += MULTlocal[1]*real(GJ_2[2]);
					Vimag[shift+2*NO+jj] += MULTlocal[1]*imag(GJ_2[2]);
					
				}
				
				// third edge
				if (IDXlocal[2] > 0){ 
				
					// we have non external edge
					
					// shift to the initial position of the edge related data (IDXLocal[2] in 3rd dimension)
					shift = (IDXlocal[2]-1)*3*NO; // shift 3*NO for each value of edge, recall: matlab indexes start at 1
					
					//first component
					Vreal[shift+jj] += MULTlocal[2]*real(GJ_3[0]);
					Vimag[shift+jj] += MULTlocal[2]*imag(GJ_3[0]);
					
					//second component (NO extra shift)
					Vreal[shift+NO+jj] += MULTlocal[2]*real(GJ_3[1]);
					Vimag[shift+NO+jj] += MULTlocal[2]*imag(GJ_3[1]);
					
					//third component (2*NO extra shift)
					Vreal[shift+2*NO+jj] += MULTlocal[2]*real(GJ_3[2]);
					Vimag[shift+2*NO+jj] += MULTlocal[2]*imag(GJ_3[2]);
					
				}
				

			} // end else if outvec > 0

			
		} // end for ii
		
	} // end for jj
	
	// ----------------------------------------------------------------------------
    // And it is done...
    //
   
}
