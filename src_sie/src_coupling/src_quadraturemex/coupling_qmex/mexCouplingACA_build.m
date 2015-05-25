% clear all
% close all
% clc

mex -output ompQuadCoil2ScatACA COMPFLAGS="$COMPFLAGS /openmp" -v -O -largeArrayDims QuadCoil2ScatACA_ompmex.cpp Coupling.cpp Coupling_ACA.cpp