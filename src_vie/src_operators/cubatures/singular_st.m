function [I_ST] = singular_st(Np_1D,ko,dx,dy,method)
%%

% GEOMETRY
r1 = [0, 0,0]';
r2 = [dx,0,0]';
r3 = [dx,dy,0]';
r4 = [0 ,dy,0]';

% Evaluate ST

if strcmp(method ,'DIRECTFN')

    [I_T1_T1] = mexdirectfn_ws_st_const(r1,r2,r3,ko,Np_1D);
    [I_T1_T2] = mexdirectfn_ws_ea_const(r1,r3,r2,r4, ko, Np_1D);
    [I_T2_T1] = mexdirectfn_ws_ea_const(r1,r3,r4,r2, ko, Np_1D);
    [I_T2_T2] = mexdirectfn_ws_st_const(r1,r3,r4,ko,Np_1D);
    
    %
    I_ST = I_T1_T1 + I_T1_T2 + I_T2_T1 + I_T2_T2;
    
elseif strcmp(method ,'DEMCEM')

    [w,z] = gauss_1d(Np_1D);
    %
    N_theta = Np_1D;
    N_psi = Np_1D;
    %
    [ w_theta , z_theta ] = gauss_1d ( N_theta );
    [ w_psi , z_psi ] = gauss_1d ( N_psi );
    
    % Evaluate ST
    
    [I_T1_T1] = mexdirect_ws_st_const(r1,r2,r3,ko,Np_1D,z,w);
    [I_T1_T2] = mexdirect_ws_ea_const(r1,r3,r2,r4, ko, N_theta, N_psi , w_theta, z_theta, w_psi, z_psi);
    [I_T2_T1] = mexdirect_ws_ea_const(r1,r3,r4,r2, ko, N_theta, N_psi , w_theta, z_theta, w_psi, z_psi);
    [I_T2_T2] = mexdirect_ws_st_const(r1,r3,r4,ko,Np_1D,z,w);
    
    I_ST = I_T1_T1 + I_T1_T2 + I_T2_T1 + I_T2_T2;
    
else
    
    fprintf('\nERROR: method should be DEMCEM or DIRECTFN')
    
end

