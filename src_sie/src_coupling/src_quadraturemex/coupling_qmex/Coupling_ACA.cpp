#include "mex.h"
#include <math.h>
#include "matrix.h"
#include "stdint.h"
#include <omp.h>
#include <complex>
#include <algorithm>

using namespace std;

void Coupling ( const double r1[3], const double r2[3], const double r3[3], const double r_obs[3], const double ko, const int Np_2D, const double Z1[], const double Z2[], const double Z3[], const double wp[], complex<double> GJ_1[3], complex<double> GJ_2[3], complex<double> GJ_3[3]  );

void Coupling_column(const mxArray* contributors, const int64_t J, const int64_t Mc,
					 const double R1[], const double R2[], const double R3[], const double RO[],
					 const double K0, const int NP_2D, const double Z1[], const double Z2[],
					 const double Z3[], const double WP[], complex<double> uJ[]) {
					 
	const int nsubs = 2;
	int64_t cell_index = 0;
	
	double r_obs[3];
	double r1[3], r2[3], r3[3];
	complex<double> GJ[3][3];
	
	for(int d = 0;d < 3;d++) {
		size_t subs[] = {J, d};
		cell_index = mxCalcSingleSubscript(contributors, nsubs, subs);
		mxArray* mx_contributors = (mxArray*) mxGetCell(contributors, cell_index);

		if(mxIsEmpty(mx_contributors))
			continue;

		double* local_contributors = (double*) mxGetData(mx_contributors);
		const int64_t num_contributors = mxGetNumberOfElements(mx_contributors);
		//mexPrintf("\nd: %d, num: %d", d, num_contributors);
		
		for(int ni=0;ni < num_contributors;ni++) {
			// get index of local contributor
			int64_t n = (int64_t) local_contributors[ni];
			double mult = (n > 0) ? 1.0 : -1.0;
			n = abs(n)-1;

			// get the coordinates of the element points
			for(int64_t r=0;r < 3;r++) {
				r1[r] = R1[3*n+r];
				r2[r] = R2[3*n+r];
				r3[r] = R3[3*n+r];
			}
			// fill up column
			#pragma omp parallel for \
			default(shared) private(r_obs,GJ)
			for(int m = 0;m < Mc;m++) {
				// get the coordinates of the observation point
				for(int r = 0;r < 3;r++)
					r_obs[r] = RO[3*m+r]; // d -> x, y or z
				
				Coupling(r1, r2, r3, r_obs, K0, NP_2D, Z1, Z2, Z3, WP, GJ[0], GJ[1], GJ[2]);
				
				for(int r = 0;r < 3;r++)
					uJ[m+r*Mc] += mult*GJ[d][r];
			}
		}
	}
}

void Coupling_row(int64_t const I, int64_t const M, int64_t const Mc, int64_t const Ne, double const* RO,
				  double const* R1, double const* R2, double const* R3, double const K0,
				  int const NP_2D, double const* Z1, double const* Z2, double const* Z3,
				  double const* WP, double const* IDX, double const* MULT, complex<double>** const Vk) {
	
	double r_obs[3];
	double r1[3], r2[3], r3[3];
	int64_t IDXlocal;
	double MULTlocal;
	complex<double> GJ[3][3];
	
	// get the coordinates of the observation point
	// constant for a row
	for(int d = 0;d < 3;d++)
		r_obs[d] = RO[3*(I % Mc)+d]; // d indicates x, y or z
	
	#pragma omp parallel for default(shared) \
	private(r1,r2,r3,GJ,IDXlocal,MULTlocal) ordered schedule(dynamic)
	for(int n = 0;n < Ne;n++) {	// Ne?
		// get the coordinates of the element points
		for(int d=0;d < 3;d++) {
			r1[d] = R1[3*n+d];
			r2[d] = R2[3*n+d];
			r3[d] = R3[3*n+d];
		}
		// call coupling function -- returns GJ
		Coupling(r1, r2, r3, r_obs, K0, NP_2D, Z1, Z2, Z3, WP, GJ[0], GJ[1], GJ[2]);
		// add computed elements to Vk
		#pragma omp critical
		for(int d = 0;d < 3;d++) {
			IDXlocal = (int64_t) (IDX[3*n+d]-1);
			MULTlocal = MULT[3*n+d];
			if(IDXlocal > -1) {
				Vk[0][IDXlocal] += conj(MULTlocal*GJ[d][0]);
				Vk[1][IDXlocal] += conj(MULTlocal*GJ[d][1]);
				Vk[2][IDXlocal] += conj(MULTlocal*GJ[d][2]);
			}
		}
	}
}

/*

*/