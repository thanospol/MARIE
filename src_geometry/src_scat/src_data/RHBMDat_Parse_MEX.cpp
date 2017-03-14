#include "mex.h"
#include <math.h>
#include "matrix.h"
#include "stdint.h"
#include <cstdio>
#ifndef __APPLE__
#include <omp.h>
#endif
#include <complex>
#include <algorithm>
#include <vector>
#include <cstring>

// Inputs: voxelfile, xfile, yfile, zfile
// Outputs: x, y, z, epsilon_r, sigma_e, rho, mu_r, sigma_m, material

using namespace std;

char seekEndOfLine(FILE* const);
void moveBack(FILE* const, int const);
void readPositionFile(char* const, int*, vector<double>&);
int readDataFile(char* const, int const, int const, int const, double* data[]);

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	char* voxelfile;
	char* xfile;
	char* yfile;
	char* zfile;

	voxelfile = mxArrayToString(prhs[0]);
	xfile = mxArrayToString(prhs[1]);
	yfile = mxArrayToString(prhs[2]);
	zfile = mxArrayToString(prhs[3]);

	mexPrintf("\n Loading \"%s\" with\n\t%s\n\t%s\n\t%s\n as coordinate files...", voxelfile, xfile, yfile, zfile);
	mexEvalString("drawnow;");

	// count number of lines in each file
	vector<double> xVec(1000,0);
	vector<double> yVec(1000,0);
	vector<double> zVec(1000,0);

	int L, M, N;
	
	readPositionFile(xfile,&L,xVec);
	readPositionFile(yfile,&M,yVec);
	readPositionFile(zfile,&N,zVec);

	plhs[0] = mxCreateUninitNumericMatrix(L, 1, mxDOUBLE_CLASS, mxREAL);
	plhs[1] = mxCreateUninitNumericMatrix(M, 1, mxDOUBLE_CLASS, mxREAL);
	plhs[2] = mxCreateUninitNumericMatrix(N, 1, mxDOUBLE_CLASS, mxREAL);

	copy(xVec.begin(), xVec.begin()+L, mxGetPr(plhs[0]));
	copy(yVec.begin(), yVec.begin()+M, mxGetPr(plhs[1]));
	copy(zVec.begin(), zVec.begin()+N, mxGetPr(plhs[2]));

	size_t dims[3] = {(size_t)L,(size_t)M,(size_t)N};
	double* lhsp[6];
	for(int k=0;k < 6;k++) {
		plhs[k+3] = mxCreateUninitNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
		lhsp[k] = mxGetPr(plhs[k+3]);
	}
	// lhs:  epsilon_r(0), sigma_e(1), rho(2), mu_r(3), sigma_m(4), material(5)
	// file: material(5), epsilon_r(0), sigma_e(1), mu_r(3), sigma_h(4), rho(2)
	double* entries[6] = {lhsp[5], lhsp[0], lhsp[1], lhsp[3], lhsp[4], lhsp[2]};
	int count = readDataFile(voxelfile,L,M,N,entries);

	mexPrintf("done.\n %d lines read\n", count);
	mexEvalString("drawnow;");
}

char seekEndOfLine(FILE* const file) {
	char c;
	do
		c = fgetc(file);
	while(c != '\n' && c != EOF);
	return c;
}

void moveBack(FILE* const file, int const n) {
	fseek(file, ftell(file)-n, SEEK_SET);
}

void readPositionFile(char* const filename, int* length, vector<double>& p) {
	FILE* const file = fopen(filename, "r");
	if(file == NULL) {
		mexErrMsgTxt("Unable to locate or open coordinate file");
	}
	int count = 0;
	double val;
	char c;
	do {
		c = fgetc(file);
		if(c == '%') {
			c = seekEndOfLine(file);
		}
		else if(c != EOF) {
			moveBack(file, 1);
			if(fscanf(file, " %lf ", &val) == 0) {
				mexErrMsgTxt("Unable to parse coordinate file properly.");
			}
			p[count] = val;
			count++;
			moveBack(file, 1);
			c = seekEndOfLine(file);
		}
	} while(c != EOF);

	(*length) = count;
	fclose(file);
}

int readDataFile(char* const filename, int const L, int const M, int const N, double* data[6]) {
	FILE* const file = fopen(filename, "r");
	if(file == NULL) {
		mexErrMsgTxt("Unable to locate or open coordinate file");
	}
	int count = 0;
	double val;
	char c;
	int x, y, z, idx;
	double ldata[6];
	do {
		c = fgetc(file);
		if(c == '%') {
			c = seekEndOfLine(file);
		}
		else if(c != EOF) {
			moveBack(file, 1);
			if(fscanf(file, " %d %d %d %lf %lf %lf %lf %f %lf ",&x,&y,&z,ldata,ldata+1,ldata+2,ldata+3,ldata+4,ldata+5) != 9) {
				mexErrMsgTxt("Unable to parse data file properly.");
			}

			idx = x + y*L + z*(L*M);
			for(int k=0;k < 6;k++)
				data[k][idx] = ldata[k];
			count++;
			moveBack(file, 1);
			c = seekEndOfLine(file);
		}
	} while(c != EOF);

	if(count != L*M*N) {
		mexWarnMsgIdAndTxt("MATLAB:RHBMDat_Parse:SizeMismatch", "Counted %d voxels when loading, but %d positions were specified", count, L*M*N);
	}

	fclose(file);
	return count;
}
