function [Zc2h] = Assembly_WCOUP(Scoord,Ccoil,Dcoil,freq) 
%%    Generate the coupling matrix for the VIE+WIE solver
% _________________________________________________________________________
%
%   Fucntion to generate the POINT MATCHING Coupling MoM to VIE
%   It is based on point matching, so quantities are point fields
%   (not voxel based fields)
%
%   Uses dyadic Green functions
%   subfunction: E_field_DGF
%
% _________________________________________________________________________
%
%%  Input
%       Scoord - coordinates of the observation points (No x 3)
%       Ccoil - coordinates of the centers of coil segments 
%       Dcoil - Cartesian components of the segments
%       freq - frequency
%
%%  Output
%       Zc2h - Tensor (No x 3 x Nd) with the contribution of each segment
%              Zc2h(1000,3,200) is z contribution of 200-th segment to 1000-th element 
%
% -------------------------------------------------------------------------
%
%   J. Fernandez Villena -- jvillena@mit.edu
%   A.G. Polimeridis -- thanos_p@mit.edu
%   Computational Prototyping Group, RLE at MIT
%
% _________________________________________________________________________

% fid = 1;

% -------------------------------------------------------------------------
%            Define EM constants
% -------------------------------------------------------------------------

co = 299792458;
lambda  = co/freq;
ko = 2*pi/lambda;

% -------------------------------------------------------------------------
% Define variables and allocate space
% -------------------------------------------------------------------------

No = size(Scoord,1); % number of observation points
Ns = size(Ccoil,1); % number of segments
Zc2h = zeros(No,3,Ns);

% -------------------------------------------------------------------------
% loop on the elements and fill the matrix
% -------------------------------------------------------------------------

% tic;

% loop on the number of elements and evaluate the effect
parfor jj = 1:Ns
    
    % D(jj,:) are the components of the current
    % S are the coordinates of DEIM points
    % C are the coordinates of the center of the segments
    Jj = Dcoil(jj,:).';
    Ei = E_field_DGF(Jj,Ccoil(jj,:),Scoord,ko);
    Zc2h(:,:,jj) = Ei;
    
end

% Timedyad = toc;

% -------------------------------------------------------------------------
%             and it is done: report 
% -------------------------------------------------------------------------


% fprintf(fid, '\n Coupling done!');
% fprintf(fid, '\n  Time in Dyad eval: %.2f', Timedyad);
% fprintf(fid, '\n ----------------------------------------------------------\n\n ');
    
