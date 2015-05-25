function [REG] = regression_svie_sparameters(UPDATE)

%%    Regression Script for solving the S-parameters of a single 8-port coil 
% _________________________________________________________________________
%
%       Script to test the SIE solver with a 8 port planar coil
%       Runs with and without an homogeneous cube
%
% _________________________________________________________________________
%
%
% -------------------------------------------------------------------------
%
%%   This function is part of MARIE
%   MARIE - Magnetic Resonance Integral Equation suite
%           Jorge Fernandez Villena   -- jvillena@mit.edu
%           Athanasios G. Polimeridis -- thanos_p@mit.edu
%           Copyright © 2014
%           RLE Computational Prototyping Group, MIT
% 
%           This software is free and open source
%           Distributed under the GNU-GPLv3 terms
%           For details see MARIE_license.txt
%
% _________________________________________________________________________


T_elapsed = tic;


% -------------------------------------------------------------------------
% Load the data from the golden folder
% -------------------------------------------------------------------------

if ispc
    load('.\golden\PLANAR_1coil_8ports_golden.mat');
else
    load('./golden/PLANAR_1coil_8ports_golden.mat');
end

% -------------------------------------------------------------------------
% Compute the Frequency response
% -------------------------------------------------------------------------

NPORT = length(SCOIL.port);

ZP_sie = zeros(NPORT,NPORT,length(fHz));
ZP_iter = zeros(NPORT,NPORT,length(fHz));


for kk = 1:length(fHz)
    
    freq = fHz(kk);
    
    % -------------------------------------------------------------------------
    % Solve the free space case
    % -------------------------------------------------------------------------
    
    titer = tic;
    
    [ZP] = SIE_Solver(SCOIL,freq);
    
    ZP_sie(:,:,kk) = ZP;
    
    fprintf(1, '\n  -- Free space solver applied');
    fprintf(1, '\n          Elapsed time  = %g [sec]', toc(titer));
    fprintf(1, '\n');

    % -------------------------------------------------------------------------
    % Solve using the iterative scheme
    % -------------------------------------------------------------------------
    
    gpu_flag = 1;
        
    titer = tic;
        
    [ZPiter,Jciter,Jbiter] = VSIE_Solver(RHBM,SCOIL,freq,tol,coup,gpu_flag);
    
    ZP_iter(:,:,kk) = ZPiter;
    
    fprintf(1, '\n  -- Iterative method applied');
    fprintf(1, '\n          Elapsed time  = %g [sec]', toc(titer));
    fprintf(1, '\n');
    
    
    
end


Z0 = 50;
SP_sie = z2s(ZP_sie, Z0);
SP_iter = z2s(ZP_iter, Z0);


% -------------------------------------------------------------------------
% validate
% -------------------------------------------------------------------------

% overall time
T_elapsed = toc(T_elapsed); 

if UPDATE
    
    % store solutions as new golden
    ZP_sie_golden = ZP_sie;
    SP_sie_golden = SP_sie;
    ZP_iter_golden = ZP_iter;
    SP_iter_golden = SP_iter;
    Jciter_golden = Jciter;
    Jbiter_golden = Jbiter;
    T_elapsed_golden = T_elapsed;
    
    if ispc
        save('.\golden\PLANAR_1coil_8ports_golden.mat', 'RHBM', 'SCOIL', 'fHz', 'tol', 'coup',...
            'ZP_sie_golden', 'SP_sie_golden', 'ZP_iter_golden', 'SP_iter_golden', 'ZP_sie_golden', 'T_elapsed_golden',...
            'Jciter_golden', 'Jbiter_golden', '-v7.3');
    else
        save('./golden/PLANAR_1coil_8ports_golden.mat', 'RHBM', 'SCOIL', 'fHz', 'tol', 'coup',...
            'ZP_sie_golden', 'SP_sie_golden', 'ZP_iter_golden', 'SP_iter_golden', 'ZP_sie_golden', 'T_elapsed_golden',...
            'Jciter_golden', 'Jbiter_golden', '-v7.3');
    end
    
    REG.TESTCASE = 'S parameters 8-port coil with homogeneous cube';
    REG.NAME{1} = 'Update';
    REG.ERROR(1) = 0;
    REG.TNEW = T_elapsed;
    REG.TOLD = T_elapsed;
    
else
    
    % Compare to golden    
    
    % TODO: check if metric makes sense
    REG.TESTCASE = 'S parameters 8-port coil with homogeneous cube';
    REG.NAME{1} = 'Z Parameters SIE';
    REG.ERROR(1) = 2*norm(ZP_sie_golden(:) - ZP_sie(:))/norm(ZP_sie_golden(:)+ZP_sie(:));
    REG.NAME{2} = 'S Parameters SIE';
    REG.ERROR(2) = 2*norm(SP_sie_golden(:) - SP_sie(:))/norm(SP_sie_golden(:) + SP_sie(:));
    REG.NAME{3} = 'Z Parameters SIE + VIE Iterative';
    REG.ERROR(3) = 2*norm(ZP_iter_golden(:) - ZP_iter(:))/norm(ZP_iter_golden(:) + ZP_iter(:));
    REG.NAME{4} = 'S Parameters SIE + VIE Iterative';
    REG.ERROR(4) = 2*norm(SP_iter_golden(:) - SP_iter(:))/norm(SP_iter_golden(:) + SP_iter(:));
    REG.NAME{5} = 'J coil SIE + VIE Iterative';
    REG.ERROR(5) = 2*norm(Jciter_golden(:) - Jciter(:))/norm(Jciter_golden(:) + Jciter(:));
    REG.NAME{6} = 'J body SIE + VIE Iterative';
    REG.ERROR(6) = 2*norm(Jbiter_golden(:) - Jbiter(:))/norm(Jbiter_golden(:) + Jbiter(:));
    REG.TNEW = T_elapsed;
    REG.TOLD = T_elapsed_golden;
 
end

