% clear all
% close all
% clc

% openMPI versions... if openmp option is not available in the compiler, this will fail 
try
    mex -output ompQuadCoil2Scat COMPFLAGS="$COMPFLAGS /openmp" -v -O -largeArrayDims QuadCoil2Scat_ompmex.cpp Coupling.cpp
    mex -output ompQuadScat2Coil COMPFLAGS="$COMPFLAGS /openmp" -v -O -largeArrayDims QuadScat2Coil_ompmex.cpp Coupling.cpp
catch me
    % use the sequential version
    fprintf(1, '\n\n\n  Warning: no openmp compatible compiler\n  Mex functions will be executed in single core\n\n');
    pause(3);
end
mex -output mexQuadCoil2Scat -v -O -largeArrayDims QuadCoil2Scat_mex.cpp Coupling.cpp
mex -output mexQuadScat2Coil -v -O -largeArrayDims QuadScat2Coil_mex.cpp Coupling.cpp