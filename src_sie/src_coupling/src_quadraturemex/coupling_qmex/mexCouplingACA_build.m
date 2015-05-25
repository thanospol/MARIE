% clear all
% close all
% clc

try
	mex -output ompQuadCoil2ScatACA COMPFLAGS="$COMPFLAGS /openmp" -v -O -largeArrayDims QuadCoil2ScatACA_ompmex.cpp Coupling.cpp Coupling_ACA.cpp
catch me
	fprintf(1, '\n\n\n  Warning: no openmp compatible compiler\n  ACA is only available on openmp-compatible clients\n\n');
	pause(3);
end