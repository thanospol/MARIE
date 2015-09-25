#include "mex.h"
#include "matrix.h"
#include "blas.h"
#include "stdint.h"
#include <cmath>
#include <omp.h>
#include <complex>
#include <algorithm>
#include <cstring>
#include <vector>

#ifndef _OPENMP
	#define omp_get_num_procs() 0
	#define omp_set_num_threads(int) 0
#endif

typedef std::complex<double> cdouble;

using namespace std;

void Coupling_column(const int64_t dtoe[], const int64_t J, const int64_t Mc, const double R1[], const double R2[], const double R3[], const double RO[], const double K0, const int NP_2D, const double Z1[], const double Z2[], const double Z3[], const double WP[], cdouble* const uJ);
void Coupling_row(int64_t const I, int64_t const M, int64_t const Mc, int64_t const Ne, double const RO[],double const R1[], double const R2[], double const R3[], double const K0,int const NP_2D, double const Z1[], double const Z2[], double const Z3[],double const WP[], int64_t const etod[], cdouble* const Vk[3]);
void update_norm(double* SF2, double* RF2, double* uvF2, int64_t M, int64_t N, int64_t k, int64_t D, double const* uk, double const* vk, double const* U, double const* V);
int compare_norm(const cdouble x, const cdouble y) {
	return norm(x) <= norm(y);
}

// pointers for BLAS
static int64_t const D0[] = {0, 1, 2, 3};
static cdouble const F[] = {-1, 0, 1, 2};
static double const* const F0 = (double* const)(F+1);

double rzdotu(int64_t const k, double* const Uu, double* const Vv) {
	#ifndef _WIN32
		return zdotu(&k, Uu, D0+1, Vv, D0+1).r;
	#else
		double UuVv[2] = {0.0,0.0};
		zgemv("N", D0+1, &k, F0+2, Uu, D0+1, Vv, D0+1, F0, UuVv, D0+1);
		return UuVv[0];
	#endif
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
// --------------------------------------
//         map IO to variables
// --------------------------------------

// check for proper datatypes
if(!mxIsInt64(prhs[6])) {
	mexErrMsgTxt("etod must be of type int64");
}
if(!mxIsInt64(prhs[14])) {
	mexErrMsgTxt("dtoe must be of type int64");
}
if(!mxIsInt64(prhs[16])) {
	mexErrMsgTxt("order must be of type int64");
}


//  point to data in vectors R 
double* const R1 = mxGetPr(prhs[0]);
double* const R2 = mxGetPr(prhs[1]);
double* const R3 = mxGetPr(prhs[2]);

// edges, observations points
int64_t const Ne = (int64_t) mxGetScalar(prhs[3]);
double* const RO = mxGetPr(prhs[4]);
int64_t const No = (int64_t) mxGetScalar(prhs[5]);

// stores sign, multiplying factor, etc.
int64_t* const etod = (int64_t*)mxGetData(prhs[6]);

int64_t const NC = (int64_t) mxGetScalar(prhs[7]);
int64_t const Nd = NC;
double const K0 = mxGetScalar(prhs[8]);
int const NP_2D = (int) mxGetScalar(prhs[9]);

//  create a pointer to w,u,v 
double* const Z1 = mxGetPr(prhs[10]);
double* const Z2 = mxGetPr(prhs[11]);
double* const Z3 = mxGetPr(prhs[12]);
double* const WP = mxGetPr(prhs[13]);

// ACA parameters
int64_t const* dtoe = (int64_t*)mxGetData(prhs[14]);
double const tol = mxGetScalar(prhs[15]);
int64_t const order = (int64_t)mxGetScalar(prhs[16]);

// --------------------------------------
//         declare local variables
// --------------------------------------

// omp
int NPROC = 1;

// constants
int64_t const M = 3*No;
int64_t const Mc = No;
int64_t const N = NC;
int64_t const D = order;

// loop variables
int64_t I = 0;
int64_t J = 0;
int64_t R = 0;

// norms
double SF2 = 0.0;
double RF2 = 0.0;
double* const uvF2 = new double[D+1];
memset(uvF2, 0, (D+1)*sizeof(double));

// components to be added to the basis
double* const uk = new double[2*M];
double* const vk = new double[2*N];
double* const vk3 = new double[3*2*N];
double* const Vk[3] = {vk3, vk3+2*N, vk3+4*N};

// cast-typed pointers for coupling calculation
cdouble* const ukc = (cdouble*)uk;
cdouble* const Vkc[3] = {(cdouble*)Vk[0], (cdouble*)Vk[1], (cdouble*)Vk[2]};

// normalization multiplier
cdouble nVkc;
double* const nVk = (double*)&nVkc;

// output arguments
vector<double> Uv;
vector<double> Vv;

double* U;
double* V;

Uv.reserve(max(2*M*min(M,N)/8, 2*M));
Vv.reserve(max(2*N*min(M,N)/8, 2*N));

// get max number of threads for omp
NPROC = omp_get_num_procs(); // get number of processors in machine
omp_set_num_threads(NPROC); // set number of threads to maximum

// loop
int64_t k = 0;
for(k = 0;k < min(M,N);k++) {
	// zero out buffers
	memset(uk,  0,   2*M*sizeof(double));
	memset(vk,  0,   2*N*sizeof(double));
	memset(vk3, 0, 3*2*N*sizeof(double));
	
	Uv.resize((k+1)*2*M);
	Vv.resize((k+1)*2*N);

	U = &Uv.front();
	V = &Vv.front();

	// calculate v_k
	Coupling_row(I, M, Mc, Ne, RO, R1, R2, R3, K0, NP_2D, Z1, Z2, Z3, WP, etod, Vkc);

	// calculate rows of R from rows of A
	for(int64_t r=0, Ir=I%Mc;r < 3;r++, Ir += Mc)
		zgemm("N", "C", &N, D0+1, &k, F0-2, V, &N, U+2*Ir, &M, F0+2, Vk[r], &N);

	// find row with largest element
	J = distance(Vkc[0], max_element(Vkc[0], Vkc[0]+3*N, compare_norm));
	R = J / N;
	J = J % N;
	
	// set new vk to row with largest element
	nVkc = 1.0/Vkc[R][J];
	zaxpy(&N, nVk, Vk[R], D0+1, vk, D0+1);

	// calculate u_k
	Coupling_column(dtoe, J, Mc, R1, R2, R3, RO, K0, NP_2D, Z1, Z2, Z3, WP, ukc);

	// calculate column of R from column of A
	zgemm("N", "C", &M, D0+1, &k, F0-2, U, &M, V+2*J, &N, F0+2, uk, &M);

   	// update norms
	update_norm(&SF2, &RF2, uvF2, M, N, k, D, uk, vk, U, V);
	
	// update basis
	zcopy(&M, uk, D0+1, U+2*k*M, D0+1);
	zcopy(&N, vk, D0+1, V+2*k*N, D0+1);
	
	// terminate or get index of new row to expand
	if(sqrt(RF2) < tol*sqrt(SF2)) {
		k += 1;
		break;
	}
	
	// find row with largest element
	ukc[I] = 0;
	I = distance(ukc, max_element(ukc, ukc+M, compare_norm));
}

// outputs
plhs[0] = mxCreateCellMatrix(1, k);
plhs[1] = mxCreateCellMatrix(1, k);

mxArray* const Uc = plhs[0];
mxArray* const Vc = plhs[1];

for(int64_t r = k-1;r >= 0;r--) {
	mxArray* Ucv = mxCreateUninitNumericMatrix(M, 1, mxDOUBLE_CLASS, mxCOMPLEX);
	mxArray* Vcv = mxCreateUninitNumericMatrix(N, 1, mxDOUBLE_CLASS, mxCOMPLEX);
	dcopy(&M, U+(2*M*r),   D0+2, mxGetPr(Ucv), D0+1);
	dcopy(&M, U+(2*M*r)+1, D0+2, mxGetPi(Ucv), D0+1);
	dcopy(&N, V+(2*N*r),   D0+2, mxGetPr(Vcv), D0+1);
	dcopy(&N, V+(2*N*r)+1, D0+2, mxGetPi(Vcv), D0+1);
	mxSetCell(Uc, r, Ucv);
	mxSetCell(Vc, r, Vcv);
	Uv.resize(r*2*M);
	Vv.resize(r*2*N);
}

}

void update_norm(double* const SF2, double* const RF2, double* const uvF2, int64_t const M, int64_t const N, int64_t const k, int64_t const D, double const* const uk, double const* const vk, double const* const U, double const* const V) {
	double uF2, vF2, UuVv;
	double* const Uu = new double[2*k];
	double* const Vv = new double[2*k];
	double* const Uu0 = new double[2*D-2];
	double* const Vv0 = new double[2*D-2];
	double UuVv0, UuVvD;
	int64_t const Dp = min(k, D-1);
	int64_t const EN0 = min(max(k-(D-1),int64_t(0)),int64_t(1));
	int64_t const Sp = EN0*(D-1);

	// blas-specific pointers
	uF2 = pow(dznrm2(&M, uk, D0+1), 2);
	vF2 = pow(dznrm2(&N, vk, D0+1), 2);

	zgemv("C", &M, &k, F0+2, U, &M, uk, D0+1, F0, Uu, D0+1);
	zgemv("C", &N, &k, F0+2, V, &N, vk, D0+1, F0, Vv, D0+1);
	zgemv("C", &M, &Sp, F0+2, U+2*(k-Dp)*M, &M, U+2*(k-Dp-1)*M, D0+1, F0, Uu0, D0+1);
	zgemv("C", &N, &Sp, F0+2, V+2*(k-Dp)*N, &N, V+2*(k-Dp-1)*N, D0+1, F0, Vv0, D0+1);

	UuVv  = 2*rzdotu(k, Uu, Vv);
	UuVvD = 2*rzdotu(Dp, Uu+2*(k-Dp), Vv+2*(k-Dp));
	UuVv0 = 2*rzdotu(Sp, Uu0, Vv0);
	
	dcopy(&D, uvF2+1, D0+1, uvF2, D0+1);
	uvF2[D] = uF2*vF2;

	(*RF2) += uvF2[D] - uvF2[0] + UuVvD - UuVv0;
	(*SF2) += uvF2[D] + UuVv;
}
