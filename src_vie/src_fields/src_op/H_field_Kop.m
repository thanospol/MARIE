function [Htot] = H_field_Kop(J,fK,Gram,Hinc,GPU_flag)
%%   Generation of H field 
% _________________________________________________________________________
%
%       Applies the K operator to a current density vector J
%       and adds the incident Magnetic field to generate
%       the total fields
%
% _________________________________________________________________________
%
%% INPUT
%   J:          Volumetric current density (LxMxNx3)
%   fK:         FFT-Circulant of K operator
%   Gram:       Gram matrix
%   freq:       frequency (in Hz) 
%   Hinc:       incident field (LxMxNx3) -- point matching!
%   GPU_flag:   applies GPU if 1, no GPU if 0
%
%
%% OUTPUT
%   Htot: total field (LxMxNx3) -- point matching!
%
%
% -------------------------------------------------------------------------
%
%   J. Fernandez Villena -- jvillena@mit.edu
%   A.G. Polimeridis -- thanos_p@mit.edu
%   Computational Prototyping Group, RLE at MIT
%
% _________________________________________________________________________


% -------------------------------------------------------------------------
% Initialization of variables
% -------------------------------------------------------------------------

if(nargin < 3 )
   fprintf(1, '\n ERROR: not enough arguments\n');
   return
end
if((nargin < 4) || isempty(Hinc))
    Hinc = 0*J; % zero incident field
end
if((nargin < 5) || isempty(GPU_flag))
    GPU_flag = 0; % no GPU
end

% Check if currents are zero for fast return
if (nnz(J) == 0)
    Htot = Hinc;
    return
end

% -------------------------------------------------------------------------
% Prepare data
% -------------------------------------------------------------------------

% domain dimensions
[L, M, N, ~] = size(J);

% positions is all the domain
idx = 1:3*L*M*N;

% ---------------------------------------------------------------------
% Call the K operator
% ---------------------------------------------------------------------

% we pass J instead of r, since r in mv_K is only use to get dimensions
[Htot] = mv_K(J, fK, J, idx, GPU_flag);
Htot = reshape(Htot,L,M,N,3);

% -------------------------------------------------------------------------
% Add the incident field to the radiated field (stored in Htot)
% -------------------------------------------------------------------------

Htot = Hinc + Htot./Gram ;

