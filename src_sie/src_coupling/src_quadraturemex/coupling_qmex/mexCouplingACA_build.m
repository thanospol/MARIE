% clear all
% close all
% clc

try
    if ispc
        mex -output ompQuadCoil2ScatACA COMPFLAGS="$COMPFLAGS /openmp" -v -O -lmwblas -largeArrayDims QuadCoil2ScatACA_ompmex.cpp Coupling.cpp Coupling_ACA.cpp;
	else
        mex -output ompQuadCoil2ScatACA COMPFLAGS='$COMPFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' CXXOPTIMFLAGS='$CXXOPTIMFLAGS -fopenmp' -lmwblas -v -O -largeArrayDims QuadCoil2ScatACA_ompmex.cpp Coupling.cpp Coupling_ACA.cpp
    end
catch me
	warning('Unable to compile ACA with openmp; make sure that your openmp libraries are installed');
	me.message
	me.stack
end