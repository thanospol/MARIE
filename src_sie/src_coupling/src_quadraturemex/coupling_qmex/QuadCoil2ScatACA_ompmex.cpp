#include "mex.h"
#include <math.h>
#include "matrix.h"
#include "stdint.h"
#include <omp.h>
#include <complex>
#include <algorithm>

using namespace std;

void Coupling_column(const mxArray* contributors, const int64_t J, const int64_t Mc, const double R1[], const double R2[], const double R3[], const double RO[], const double K0, const int NP_2D, const double Z1[], const double Z2[], const double Z3[], const double WP[], complex<double> uJ[]);
void Coupling_row(int64_t const I, int64_t const M, int64_t const Mc, int64_t const Ne, double const* RO,double const* R1, double const* R2, double const* R3, double const K0,int const NP_2D, double const* Z1, double const* Z2, double const* Z3,double const* WP, double const* IDX, double const* MULT, complex<double>** const Vk);
void update_norm(double* const SF2, double* const uvF2, int const M, int const N, int const k, complex<double> const* uk, complex<double> const* vk, double** const* U, double** const* V);
int compare_norm(const complex<double> x, const complex<double> y) {
	return norm(x) <= norm(y);
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
// --------------------------------------
//         map IO to variables
// --------------------------------------

//  point to data in vectors R 
double* const R1 = mxGetPr(prhs[0]);
double* const R2 = mxGetPr(prhs[1]);
double* const R3 = mxGetPr(prhs[2]);

// edges, observations points
int64_t const Ne = (int64_t) mxGetScalar(prhs[3]);
double* const RO = mxGetPr(prhs[4]);
int64_t const No = (int64_t) mxGetScalar(prhs[5]);

// stores sign, multiplying factor, etc.
double* const IDXdouble = mxGetPr(prhs[6]);
double* const MULT = mxGetPr(prhs[7]);

int64_t const NC = (int64_t) mxGetScalar(prhs[8]);
int64_t const Nd = NC;
double const K0 = mxGetScalar(prhs[9]);
int const NP_2D = (int) mxGetScalar(prhs[10]);

//  create a pointer to w,u,v 
double* const Z1 = mxGetPr(prhs[11]);
double* const Z2 = mxGetPr(prhs[12]);
double* const Z3 = mxGetPr(prhs[13]);
double* const WP = mxGetPr(prhs[14]);

// ACA parameters
const mxArray* contributors = prhs[15];
double const tol = mxGetScalar(prhs[16]);
int order = (int) mxGetScalar(prhs[17]);

// --------------------------------------
//         declare local variables
// --------------------------------------

// omp
int NPROC = 1;

// constants
const int64_t M = 3*No;
const int64_t Mc = No;
const int64_t N = NC;
const int64_t D = order;

// loop variables
int64_t I = 0;
int64_t J = 0;

// intermediate storage
double SF2 = 0.0;
double uvF2 = 0.0;

// local storage of row(s) & column
complex<double>* uk = new complex<double>[M];
complex<double>* vk = new complex<double>[N];
complex<double>** Vk = new complex<double>*[3];
for(int r = 0;r < 3;r++)
	Vk[r] = new complex<double>[N];

// output arguments
double*** const U = new double**[min(M,N)]; // 0 real
double*** const V = new double**[min(M,N)]; // 1 imag
for(int r = 0;r < min(M,N);r++) {
	U[r] = new double*[2];
	V[r] = new double*[2];
}

// get max number of threads for omp
NPROC = omp_get_num_procs(); // get number of processors in machine
omp_set_num_threads(NPROC); // set number of threads to maximum

// loop
int k;
for(k = 0;k < min(M,N);k++) {
	// --------------------------------------
	//           zero out buffers
	// --------------------------------------
	for(int m=0;m < M;m++)
		uk[m] = 0;
	for(int n=0;n < N;n++) {
		vk[n] = 0;
		for(int r=0;r < 3;r++)
			Vk[r][n] = 0;
	}

    // --------------------------------------
	//         calculate v_k
	// --------------------------------------
	
	// calculate rows of A corresponding to I
	Coupling_row(I, M, Mc, Ne, RO, R1, R2, R3, K0, NP_2D, Z1, Z2, Z3, WP, IDXdouble, MULT, Vk);

	// calculate rows of R from rows of A
	#pragma omp parallel for
	for(int64_t n = 0;n < N;n++) {
		for(int64_t r = 0;r < 3;r++) {
			int64_t Ir = (I%Mc)+r*Mc;
			complex<double> vP = 0.0;
			for(int64_t l = 0;l < k;l++)
				vP += complex<double>(V[l][0][n],V[l][1][n])*complex<double>(U[l][0][Ir],-U[l][1][Ir]);
			Vk[r][n] -= vP;
		}
	}

	// find row with largest element
	int R = 0;
	double vM = 0.0;
	for(int r = 0;r < 3;r++) {
		complex<double>* vPtr = max_element(Vk[r], Vk[r]+N, compare_norm);
		if(norm(*vPtr) > vM) {
			J = distance(Vk[r], vPtr);
			R = r;
			vM = norm(*vPtr);
		}
	}
	
	// set new vk to row with largest element
	memcpy(vk, Vk[R], N*sizeof(complex<double>));
	#pragma omp parallel for
	for(int n = 0;n < N;n++) {
		vk[n] /= Vk[R][J];
	}

    // --------------------------------------
	//         calculate u_k
	// --------------------------------------
	Coupling_column(contributors, J, Mc, R1, R2, R3, RO, K0, NP_2D, Z1, Z2, Z3, WP, uk);

	// calculate column of R from column of A
	#pragma omp parallel for
	for(int m = 0;m < M;m++) {
		complex<double> uP = 0.0;
		for(int r = 0;r < k;r++)
			uP += complex<double>(U[r][0][m], U[r][1][m])*complex<double>(V[r][0][J], -V[r][1][J]); // U[r][m]*conj(V[r][J]);
		uk[m] -= uP;
	}

    // --------------------------------------
	//        update norms and basis
	// --------------------------------------
	update_norm(&SF2, &uvF2, M, N, k, uk, vk, U, V);
	
	U[k] = new double*[2];
	V[k] = new double*[2];
	U[k][0] = (double*)mxMalloc(M*sizeof(double)); // new double[M];
	U[k][1] = (double*)mxMalloc(M*sizeof(double)); // new double[M];
	V[k][0] = (double*)mxMalloc(N*sizeof(double)); // new double[N];
	V[k][1] = (double*)mxMalloc(N*sizeof(double)); // new double[N];
	
	for(int m = 0;m < M;m++) {
		U[k][0][m] = real(uk[m]);
		U[k][1][m] = imag(uk[m]);
	}
	for(int n = 0;n < N;n++) {
		V[k][0][n] = real(vk[n]);
		V[k][1][n] = imag(vk[n]);
	}
	// memcpy(U[k], uk, M*sizeof(complex<double>));
	// memcpy(V[k], vk, N*sizeof(complex<double>));
	
	// ---------------------------------------------
	//  terminate or get index of new row to expand
	// ---------------------------------------------
	if(sqrt(uvF2) < tol*sqrt(SF2)) {
		k += 1;
		break;
	}
	
	// find row with largest element
	uk[I] = 0;
	I = distance(uk, max_element(uk, uk+M, compare_norm));
}

// outputs
plhs[0] = mxCreateCellMatrix(1, k);
plhs[1] = mxCreateCellMatrix(1, k);

mxArray* const Uc = plhs[0];
mxArray* const Vc = plhs[1];

//#pragma omp parallel for
const mwSize Udims[2] = {M, 1};
const mwSize Vdims[2] = {N, 1};

for(int r = 0;r < k;r++) {
	mxArray* Ucv = mxCreateNumericArray(2, Udims, mxDOUBLE_CLASS, mxCOMPLEX);
	mxArray* Vcv = mxCreateNumericArray(2, Vdims, mxDOUBLE_CLASS, mxCOMPLEX);
	mxSetPr(Ucv, U[r][0]);
	mxSetPi(Ucv, U[r][1]);
	mxSetPr(Vcv, V[r][0]);
	mxSetPi(Vcv, V[r][1]);
	mxSetCell(Uc, r, Ucv);
	mxSetCell(Vc, r, Vcv);
}

}

void update_norm(double* const SF2, double* const uvF2, int const M, int const N, int const k, complex<double> const* uk, complex<double> const* vk, double** const* U, double** const* V) {
	double uF2 = 0.0;
	double vF2 = 0.0;
	
	double UVxF2 = 0.0;
	complex<double> Uu = 0.0;
	complex<double> Vv = 0.0;
	
	for(int m=0;m < M;m++)
		uF2 += norm(uk[m]);

	for(int n=0;n < N;n++)
		vF2 += norm(vk[n]);

	#pragma omp parallel for reduction(+:UVxF2) \
	 default(shared) private(Uu,Vv)
	for(int r=0;r < k;r++) {
		Uu = 0.0;
		for(int m = 0;m < M;m++)
			Uu += complex<double>(U[r][0][m],-U[r][1][m])*uk[m]; // conj(U[r][m])*uk[m];
		Vv = 0.0;
		for(int n = 0;n < N;n++)
			Vv += complex<double>(V[r][0][n],-V[r][1][n])*vk[n]; // conj(V[r][n])*vk[n];
		UVxF2 += 2*real(Uu*Vv);
	}

	(*uvF2) = uF2*vF2;
	(*SF2) += (*uvF2) + UVxF2;
}