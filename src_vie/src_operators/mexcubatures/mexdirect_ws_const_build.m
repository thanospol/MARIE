
% Nop
mex -output mexcubature_nop -v -O cubature_nop_mex.cpp cubature_nop.cpp
% Kop
mex -output mexcubature_kop -v -O cubature_kop_mex.cpp cubature_kop.cpp 
% DEMCEM
mex -output mexdirect_ws_st_const -largeArrayDims -v -O direct_ws_st_const_mex.cpp direct_ws_st_const.cpp
mex -output mexdirect_ws_ea_const -largeArrayDims -v -O direct_ws_ea_const_mex.cpp direct_ws_ea_const.cpp 
mex -output mexdirect_ws_va_const -largeArrayDims -v -O direct_ws_va_const_mex.cpp direct_ws_va_const.cpp 
% DIRECTFN
mex -output mexdirectfn_ws_st_const -largeArrayDims -v -O directfn_ws_st_const_mex.cpp directfn_ws_st_const.cpp
mex -output mexdirectfn_ws_ea_const -largeArrayDims -v -O directfn_ws_ea_const_mex.cpp directfn_ws_ea_const.cpp 
mex -output mexdirectfn_ws_va_const -largeArrayDims -v -O directfn_ws_va_const_mex.cpp directfn_ws_va_const.cpp 