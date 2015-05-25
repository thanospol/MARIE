function [JOut] = mv_A(JIn, fN, fK, Mer, Mce, Mmr, Mcm, Gram, ce, cm, transp_flag, idx, GPU_flag)
%%    M-V product of the complete operator via FFT
% _________________________________________________________________________
%
%       Applies the N and K operators to electric and magnetic currents
%       The currents are stored in a vector, first electric, then magnetic
%       Only the non-air components (idx) are stored
%
% _________________________________________________________________________
%
%% INPUT
%   JIn:            Current
%   fN:             FFT-Circulant of N operator
%   fK:             FFT-Circulant of K operator
%   Mer:            Relative electric permittivity
%   Mce:            Electric susceptibility 
%   Mmr:            Relative magnetic permeability
%   Mcm:            Magnetic susceptibility 
%   Gram:           Gram matrix 
%   ce:             constant jwe
%   cm:             constant jwm 
%   transp_flag:    flag for choosing the m-v product
%                   'transp':   A'* b
%                   'notransp': A * b
%   idx:            index with local coordinates (non-air voxels)
%   GPU_flag:       applies GPU if 1, no GPU if 0
%
%% OUTPUT
%   JOut:           Current
%
%
% -------------------------------------------------------------------------
%
%   A.G. Polimeridis -- thanos_p@mit.edu
%   J. Fernandez Villena -- jvillena@mit.edu
%   Computational Prototyping Group, RLE at MIT
%
% _________________________________________________________________________


% -------------------------------------------------------------------------
% Prepare data
% -------------------------------------------------------------------------

% number of non-air components
nS = length(idx);

Je = JIn(1:nS,1); % get the electric currents 
Jm = JIn((nS+1):end,1); % get the magnetic currents

% -------------------------------------------------------------------------
% Apply operations depending on transp_flag
% -------------------------------------------------------------------------

if strcmp(transp_flag,'transp')    
    
    [JeOut] = mv_AN(Je, fN, Mer, Mce, Gram, 'transp', idx, GPU_flag) - mv_AK(Je, fK, Mcm, cm, 'transp', idx, GPU_flag);
    [JmOut] = mv_AK(Jm, fK, Mce, ce, 'transp', idx, GPU_flag)        + mv_AN(Jm, fN, Mmr, Mcm, Gram, 'transp', idx, GPU_flag);
    
else
    
    [JeOut] =  mv_AN(Je, fN, Mer, Mce, Gram, 'notransp', idx, GPU_flag) + mv_AK(Jm, fK, Mce, ce, 'notransp', idx, GPU_flag);
    [JmOut] = -mv_AK(Je, fK, Mcm, cm, 'notransp', idx, GPU_flag)        + mv_AN(Jm, fN, Mmr, Mcm, Gram, 'notransp', idx, GPU_flag);
    
end

% -------------------------------------------------------------------------
% Return current vector with electric and magnetic currents
% -------------------------------------------------------------------------

JOut = [JeOut; JmOut];