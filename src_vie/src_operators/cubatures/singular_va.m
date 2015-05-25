function [I_VAc,I_VAo] = singular_va(Np_1D,ko,dx,dy,dz,method)
%%

% GEOMETRY
r1 = [0, 0,0]';
r2 = [dx,0,0]';
r3 = [dx,dy,0]';
r4 = [0 ,dy,0]';
r5 = [2*dx,dy,0]';
r6 = [2*dx ,2*dy,0]';
r7 = [dx ,2*dy,0]';
%
r5o = [dx,dy,dz]';
r6o = [dx ,2*dy,dz]';
r7o = [dx ,2*dy,0]';

%
if strcmp(method ,'DIRECTFN')
    
    [I_T1_T1c] = mexdirectfn_ws_va_const(r3,r1,r2,r5,r6, ko, Np_1D);
    [I_T1_T2c] = mexdirectfn_ws_va_const(r3,r1,r2,r6,r7, ko, Np_1D);
    [I_T2_T1c] = mexdirectfn_ws_va_const(r3,r4,r1,r5,r6, ko, Np_1D);
    [I_T2_T2c] = mexdirectfn_ws_va_const(r3,r4,r1,r6,r7, ko, Np_1D);
    
    %
    I_VAc = I_T1_T1c + I_T1_T2c + I_T2_T1c + I_T2_T2c;
    % Evaluate EAo
    [I_T1_T1o] = mexdirectfn_ws_va_const(r3,r1,r2,r5o,r6o, ko, Np_1D);
    [I_T1_T2o] = mexdirectfn_ws_va_const(r3,r1,r2,r6o,r7o, ko, Np_1D);
    [I_T2_T1o] = mexdirectfn_ws_va_const(r3,r4,r1,r5o,r6o, ko, Np_1D);
    [I_T2_T2o] = mexdirectfn_ws_va_const(r3,r4,r1,r6o,r7o, ko, Np_1D);
    
    %
    I_VAo = I_T1_T1o + I_T1_T2o + I_T2_T1o + I_T2_T2o;

elseif strcmp(method ,'DEMCEM')
    
    N_theta_p   = Np_1D;
    N_theta_q   = Np_1D;
    N_psi       = Np_1D;
    %
    [ w_theta_p,z_theta_p ] = gauss_1d ( N_theta_p );
    [ w_theta_q,z_theta_q ] = gauss_1d ( N_theta_q );
    [ w_psi,z_psi ]         = gauss_1d ( N_psi );
    % Evaluate EAc
    [I_T1_T1c] = mexdirect_ws_va_const(r3,r1,r2,r5,r6, ko, N_theta_p, N_theta_q, N_psi, w_theta_p, z_theta_p, w_theta_q, z_theta_q, w_psi, z_psi);
    [I_T1_T2c] = mexdirect_ws_va_const(r3,r1,r2,r6,r7, ko, N_theta_p, N_theta_q, N_psi, w_theta_p, z_theta_p, w_theta_q, z_theta_q, w_psi, z_psi);
    [I_T2_T1c] = mexdirect_ws_va_const(r3,r4,r1,r5,r6, ko, N_theta_p, N_theta_q, N_psi, w_theta_p, z_theta_p, w_theta_q, z_theta_q, w_psi, z_psi);
    [I_T2_T2c] = mexdirect_ws_va_const(r3,r4,r1,r6,r7, ko, N_theta_p, N_theta_q, N_psi, w_theta_p, z_theta_p, w_theta_q, z_theta_q, w_psi, z_psi);
    
    %
    I_VAc = I_T1_T1c + I_T1_T2c + I_T2_T1c + I_T2_T2c;
    % Evaluate EAo
    [I_T1_T1o] = mexdirect_ws_va_const(r3,r1,r2,r5o,r6o, ko, N_theta_p, N_theta_q, N_psi, w_theta_p, z_theta_p, w_theta_q, z_theta_q, w_psi, z_psi);
    [I_T1_T2o] = mexdirect_ws_va_const(r3,r1,r2,r6o,r7o, ko, N_theta_p, N_theta_q, N_psi, w_theta_p, z_theta_p, w_theta_q, z_theta_q, w_psi, z_psi);
    [I_T2_T1o] = mexdirect_ws_va_const(r3,r4,r1,r5o,r6o, ko, N_theta_p, N_theta_q, N_psi, w_theta_p, z_theta_p, w_theta_q, z_theta_q, w_psi, z_psi);
    [I_T2_T2o] = mexdirect_ws_va_const(r3,r4,r1,r6o,r7o, ko, N_theta_p, N_theta_q, N_psi, w_theta_p, z_theta_p, w_theta_q, z_theta_q, w_psi, z_psi);
    
    %
    I_VAo = I_T1_T1o + I_T1_T2o + I_T2_T1o + I_T2_T2o;
    
else
    
    fprintf('\nERROR: method should be DEMCEM or DIRECTFN')
    
end