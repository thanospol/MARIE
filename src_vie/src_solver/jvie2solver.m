function [Vsol] = jvie2solver(Vrhs,fVIE,fVIEGPU,tol,Vini)
%%    Solver Routine for the Electric JVIE formulation used with SCOILS
% _________________________________________________________________________
%
%       Solves the VIE problem
%       BICGSTAB no precond, JVIE 2 formulation
%       choses best option (GPU, FFTW, standard)
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



% fast return if input is zero
if (nnz(Vrhs) == 0)
    Vsol = Vrhs;
    return
end  

% -------------------------------------------------------------------------
% set some variables
% -------------------------------------------------------------------------

fid = 1;

max_it = 1000;

% -------------------------------------------------------------------------
% Call the solver
% -------------------------------------------------------------------------
        
% try GPU or mexed functions... if there are not, go for BiCGStab
ta = tic;
Vsol = [];
if ~isempty(fVIEGPU)
    try
        % GPU based BICGStab
        [Vsol, ~, relres, ~, resvec] = bicgstab(@(J)fVIEGPU(J), Vrhs, tol, max_it, [], [], Vini);
        solvert = toc(ta);
        fprintf(fid,'\n GPU BICGSTAB: Time  = %.2f [sec] for %d iterations, relative residue %g ' ,solvert,length(resvec),relres);
    catch
        Vsol = [];
    end
end
    
if isempty(Vsol) % not solved
    % try standard BICGSTAB
    [Vsol, ~, relres, ~, resvec] = bicgstab(@(J)fVIE(J), Vrhs, tol, max_it, [], [], Vini);
    solvert = toc(ta);
    fprintf(fid,'\n CPU BICGSTAB: Time  = %.2f [sec] for %d iterations, relative residue %g ' ,solvert,length(resvec),relres);
end

% -------------------------------------------------------------------------
% and done
% -------------------------------------------------------------------------


