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
	//		J - NCx3 complex numbers with current components at each observation point
	//
	//
	//  Returns:
	//      V - NC complex with each component of the field in each of the NC edges
	//
	//
	//	NOTE: 	all I/O are treated as vectors, and ordering in column major format
	//			for the coordinates, index and mult, the ordering is for ii point, x = R[3*ii], y = R[3*ii+1], z = R[3*ii+2]
	//
	// 	MATLAB CALL:
	//			[V] = mexQuadScat2Coil(R1,R2,R3,NE,RO,NO,IDX,MULT,NC,K0,NP_2D,Z1,Z2,Z3,WP,J)
	//

	// ----------------------------------------------------------------------------
	// some initial declarations
	

	if (nrhs != 16){
        mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumInputs", "QuadCoil2Scat_mex - 16 inputs required!");
	}
	
	// ----------------------------------------------------------------------------
	// declarations
	
	// aux vars
	int NPROC = 1;
	int64_t ii = 0, jj = 0, shift = 0;
	double auxreal = 0.0, auximag = 0.0;
	
	//  point to data in vectors R 
	double* R1 = mxGetPr(prhs[0]);
	double* R2 = mxGetPr(prhs[1]);
	double* R3 = mxGetPr(prhs[2]);
	int64_t NE = (int64_t) mxGetScalar(prhs[3]);
	double* RO = mxGetPr(prhs[4]);
	int64_t NO = (int64_t) mxGetScalar(prhs[5]);
	
	// pointer to the indexes
	double* IDXdouble = mxGetPr(prhs[6]); // get complex index vector!!!
	
	// pointer to the multipliers
	double* MULT = mxGetPr(prhs[7]); // get complex index vector!!!
	
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
	
	// pointer to the body currents
	double * Jreal = mxGetPr(prhs[15]);
    double * Jimag = mxGetPi(prhs[15]);

		
  	// this is for the computation of each coupling entry
    complex<double> GJ_1[3], GJ_2[3], GJ_3[3];
	double r1[3], r2[3], r3[3], r_obs[3];
	int64_t IDXlocal[3]; // index local
	double MULTlocal[3]; // local multiplier

	// for the output matrix V
	double *Vreal, *Vimag;
	
	// allocate space for the NC complex vector
	// allocate MATLAB size for the outputs, and assign pointer
    plhs[0] = mxCreateDoubleMatrix(NC, 1, mxCOMPLEX); // this initializes to 0
    Vreal = mxGetPr(plhs[0]);
    Vimag = mxGetPi(plhs[0]);


    // ----------------------------------------------------------------------------
    // We are finally good to do some computations...
    //
	
    NPROC = omp_get_num_procs(); // get number of processors in machine
	omp_set_num_threads(NPROC); // set number of threads to maximum
	
    // make sure that if the imaginary part of J is zero, this does not break
    if (Jimag == NULL){
        Jimag = (double*) calloc (3*NO,sizeof(double)); // just zero array
    }
    
    #pragma omp parallel for \
	default(shared) private(r_obs,r1,r2,r3,GJ_1,GJ_2,GJ_3,IDXlocal,MULTlocal,shift,jj,ii,auxreal,auximag)
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
	
		// get the index of the different edges
		IDXlocal[0] = (int64_t) IDXdouble[3*ii]; // index of first edge
		IDXlocal[1] = (int64_t) IDXdouble[3*ii+1]; // index of second edge
		IDXlocal[2] = (int64_t) IDXdouble[3*ii+2]; // index of third edge
		
		// get the positive or negative multiplier
		MULTlocal[0] = MULT[3*ii]; // multiplier of first edge
		MULTlocal[1] = MULT[3*ii+1]; // multiplier of second edge
		MULTlocal[2] = MULT[3*ii+2]; // multiplier of third edge
		
		
		for(jj = 0; jj < NO; jj++){ // loop on the number of observation points

			// get the coordinates of the observation point
			r_obs[0] = RO[3*jj]; // x coordinate
			r_obs[1] = RO[3*jj+1]; // y coordinate
			r_obs[2] = RO[3*jj+2]; // z coordinate		
			
			
			// Call the C-subroutine for each case
			Coupling ( r1, r2, r3, r_obs, K0, NP_2D, Z1, Z2, Z3, WP, GJ_1, GJ_2, GJ_3  );
			
			
			// apply the product according to the output case
				
			// in matlab:
			// Vout(idx1,1) = Vout(idx1,1) + mult(1)*(Gcoup(1,1)*Vin(jj,1) + Gcoup(1,2)*Vin(jj,2) + Gcoup(1,3)*Vin(jj,3)); % first edge
			

			// first edge
			if (IDXlocal[0] > 0){ 
			
				// we have non external edge
				shift = IDXlocal[0] - 1; // recall matlab indexes start in 1
				
				// first component
				auxreal = Jreal[jj]*real(GJ_1[0]) - Jimag[jj]*imag(GJ_1[0]);
				auximag = Jreal[jj]*imag(GJ_1[0]) + Jimag[jj]*real(GJ_1[0]);
				
				// second component
				auxreal += Jreal[NO+jj]*real(GJ_1[1]) - Jimag[NO+jj]*imag(GJ_1[1]);
				auximag += Jreal[NO+jj]*imag(GJ_1[1]) + Jimag[NO+jj]*real(GJ_1[1]);
				
				// third component
				auxreal += Jreal[2*NO+jj]*real(GJ_1[2]) - Jimag[2*NO+jj]*imag(GJ_1[2]);
				auximag += Jreal[2*NO+jj]*imag(GJ_1[2]) + Jimag[2*NO+jj]*real(GJ_1[2]);
				
				// weight and assign contribution to edge
                #pragma omp atomic
				Vreal[shift] += MULTlocal[0]*auxreal;
                #pragma omp atomic
				Vimag[shift] += MULTlocal[0]*auximag;
				
			}
						
			// second edge
			if (IDXlocal[1] > 0){ 
			
				// we have non external edge
				shift = IDXlocal[1] - 1; // recall matlab indexes start in 1
				
				// first component
				auxreal = Jreal[jj]*real(GJ_2[0]) - Jimag[jj]*imag(GJ_2[0]);
				auximag = Jreal[jj]*imag(GJ_2[0]) + Jimag[jj]*real(GJ_2[0]);
				
				// second component
				auxreal += Jreal[NO+jj]*real(GJ_2[1]) - Jimag[NO+jj]*imag(GJ_2[1]);
				auximag += Jreal[NO+jj]*imag(GJ_2[1]) + Jimag[NO+jj]*real(GJ_2[1]);
				
				// third component
				auxreal += Jreal[2*NO+jj]*real(GJ_2[2]) - Jimag[2*NO+jj]*imag(GJ_2[2]);
				auximag += Jreal[2*NO+jj]*imag(GJ_2[2]) + Jimag[2*NO+jj]*real(GJ_2[2]);
				
				// weight and assign contribution to edge
                #pragma omp atomic
				Vreal[shift] += MULTlocal[1]*auxreal;
                #pragma omp atomic
				Vimag[shift] += MULTlocal[1]*auximag;
				
			}
			
			// third edge
			if (IDXlocal[2] > 0){ 
			
				// we have non external edge
				shift = IDXlocal[2] - 1; // recall matlab indexes start in 1
				
				// first component
				auxreal = Jreal[jj]*real(GJ_3[0]) - Jimag[jj]*imag(GJ_3[0]);
				auximag = Jreal[jj]*imag(GJ_3[0]) + Jimag[jj]*real(GJ_3[0]);
				
				// second component
				auxreal += Jreal[NO+jj]*real(GJ_3[1]) - Jimag[NO+jj]*imag(GJ_3[1]);
				auximag += Jreal[NO+jj]*imag(GJ_3[1]) + Jimag[NO+jj]*real(GJ_3[1]);
				
				// third component
				auxreal += Jreal[2*NO+jj]*real(GJ_3[2]) - Jimag[2*NO+jj]*imag(GJ_3[2]);
				auximag += Jreal[2*NO+jj]*imag(GJ_3[2]) + Jimag[2*NO+jj]*real(GJ_3[2]);
				
				// weight and assign contribution to edge
                #pragma omp atomic
				Vreal[shift] += MULTlocal[2]*auxreal;
                #pragma omp atomic
				Vimag[shift] += MULTlocal[2]*auximag;
				
			}			

			
		} // end for jj
		
	} // end for ii
	
	// ----------------------------------------------------------------------------
    // And it is done...
    //
   
}
