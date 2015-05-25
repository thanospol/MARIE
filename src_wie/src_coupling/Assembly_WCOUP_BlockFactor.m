function [ZbcL,ZbcR] = Assembly_WCOUP_BlockFactor(Scoord,Ccoil,Dcoil,freq,tol,blockSize) 
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

% -------------------------------------------------------------------------
% loop on the elements and fill the matrix
% -------------------------------------------------------------------------

% tic;

ZbcL = [];
ZbcR = [];
svTresh = [];

iniElem = 0;
blockSize = min(blockSize,Ns-iniElem);

while (blockSize > 0)

    Zblock = zeros(No,3,blockSize);
    parfor jj = 1:blockSize

        % D(jj,:) are the components of the current
        % S are the coordinates of DEIM points
        % C are the coordinates of the center of the segments
        Jj = Dcoil(jj,:).';
        Zblock(:,:,jj) = E_field_DGF(Jj,Ccoil(jj,:),Scoord,ko);
    
    end
    Zblock = reshape(Zblock,3*No,blockSize);
    
    
    % backorth to previuos block
    if ~isempty(ZbcL)
        alpha = ZbcL'*Zblock;
        Zblock = Zblock - ZbcL*alpha;
    end
    
    % apply SVD on remaining block
    [Ublock, Sblock, Vblock] = svd(Zblock,'econ');
    clear Zblock;
    
    % define treshold for truncation
    if isempty(svTresh)
        svTresh = tol*Sblock(1,1);
    end
    % find the breaking point
    for kk = 1:length(Sblock)
        if (Sblock(kk,kk) < svTresh)
            break
        end
    end
    % truncate
    Ublock = Ublock(:,1:kk);
    Sblock = Sblock(1:kk,1:kk);
    Vblock = Vblock(:,1:kk);
    Vblock = Sblock*Vblock';
    clear Sblock;
    
    % store the new blocks in corresponding positions
    ZbcL = [ZbcL Ublock]; % append columns
    clear Ublock;
    
    % store the right factor in the corresponding positions
    % note that positions are dictated by index
    % and several elements may have the same index
    % so we have to take overlaps into account
    if ~isempty(ZbcR)
        [prevrows, prevcols] = size(ZbcR);
        newrows = size(ZbcL,2);
        newcols = prevcols + size(alpha,2);
        ZbcR(newrows,newcols) = 0; % increases size with zeros
        ZbcR(1:prevrows,prevcols+1:newcols) = alpha; % append colummns
        ZbcR(prevrows+1:newrows,prevcols+1:newcols) = Vblock; % append rows        
    else
        ZbcR = Vblock; % append rows
    end
    clear alpha; clear Vblock;

    % restart counters and arrays
    iniElem = iniElem + blockSize;
    blockSize = min(blockSize,Ns-iniElem);
    
end


% Timedyad = toc;

% -------------------------------------------------------------------------
%             and it is done: report 
% -------------------------------------------------------------------------


% fprintf(fid, '\n Coupling done!');
% fprintf(fid, '\n  Time in Dyad eval: %.2f', Timedyad);
% fprintf(fid, '\n ----------------------------------------------------------\n\n ');
    
