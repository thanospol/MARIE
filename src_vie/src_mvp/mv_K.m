function [JOut] = mv_K(JIn0, fK, r, idx, GPU_flag)
%%   M-V product of operator K via FFT
% _________________________________________________________________________
%
%       Applies the K operator to a current vector
%
% _________________________________________________________________________
%
%% INPUT
%   JIn:            Current 
%   fK:             FFT-Circulant of K operator
%   r:              domain (LxMxNx3) 
%   idx:            index with local coordinates (non-air voxels)
%   GPU_flag:       applies GPU if 1, no GPU if 0
%
%
%% OUTPUT
%   JOut:           Current
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
% Prepare data
% -------------------------------------------------------------------------

% fft dimensions
[LfK, MfK, NfK, ~] = size(fK);

% domain dimensions
[L, M, N, ~] = size(r);

% Check GPU use
if GPU_flag
    % allocate space
    JOut = gpuArray.zeros(L, M, N, 3);
    JIn = gpuArray.zeros(L, M, N, 3);    
    % send to gpu, and translate from local (idx) to global (L,M,N) coordinates
    JIn(idx) = gpuArray(JIn0);
else
    % allocate space
    JIn = zeros(L, M, N, 3);
    JOut = zeros(L, M, N, 3);    
    % translate from local (idx) to global (L,M,N) coordinates
    JIn(idx) = JIn0(:);
end

% ---------------------------------------------------------------------
% apply fft and mv-op for each of the components of JIn
% ---------------------------------------------------------------------
   
% apply fft and mv-op
fJ = fftn(JIn(:,:,:,1),[LfK, MfK, NfK]);
Jout3 =  fK(:,:,:,2) .* fJ;
Jout2 = -fK(:,:,:,3) .* fJ;
% Jout1 = 0 .* fJ;

fJ = fftn(JIn(:,:,:,2),[LfK, MfK, NfK]);
Jout1 = fK(:,:,:,3) .* fJ;
% Jout2 = Jout2 + 0 .* fJ;
Jout3 = Jout3 - fK(:,:,:,1) .* fJ;

fJ = fftn(JIn(:,:,:,3),[LfK, MfK, NfK]);
% Jout3 = Jout3 + 0 .* fJ;
Jout2 = Jout2 + fK(:,:,:,1) .* fJ;
Jout1 = Jout1 - fK(:,:,:,2) .* fJ;

% apply ifft
Jout1 = ifftn(Jout1);
JOut(:,:,:,1) = Jout1(1:L,1:M,1:N);
Jout2 = ifftn(Jout2);
JOut(:,:,:,2) = Jout2(1:L,1:M,1:N);
Jout3 = ifftn(Jout3);
JOut(:,:,:,3) = Jout3(1:L,1:M,1:N);
    
% -------------------------------------------------------------------------
% Return local coordinates related to material positions
% -------------------------------------------------------------------------
if GPU_flag
    % get from GPU
    JOut = gather(JOut(idx));
    % clear gpu data
    clear JIn; clear Jout1; clear Jout2; clear Jout3; clear fJ;
else
    JOut = JOut(idx);
end
