function [JOut] = mv_AK(JIn0, fK, Mc, c, transp_flag, idx, GPU_flag)
%%   M-V product of operator A_K via FFT
% _________________________________________________________________________
%
%       Applies the AK operator to a current vector
%
% _________________________________________________________________________
%
%% INPUT
%   JIn:            Current 
%   fK:             FFT-Circulant of K operator
%   Mc:             Susceptibility 
%   c:              constant (jwe, jwm) 
%   transp_flag:    flag for choosing the m-v product
%                   'transp':   AK' * JIn0
%                   'notransp': AK * JIn0
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
%   A.G. Polimeridis -- thanos_p@mit.edu
%   J. Fernandez Villena -- jvillena@mit.edu
%   Computational Prototyping Group, RLE at MIT
%
% _________________________________________________________________________


% -------------------------------------------------------------------------
% Prepare data
% -------------------------------------------------------------------------

% fft dimensions
[LfK, MfK, NfK, ~] = size(fK);

% domain dimensions
[L, M, N] = size(Mc);

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

% -------------------------------------------------------------------------
% Apply operations depending on transp_flag
% -------------------------------------------------------------------------

if strcmp(transp_flag,'transp')       % JOut = fK' * JIn0
    
    % ---------------------------------------------------------------------
    % apply fft and mv-op for each of the components of JIn
    % ---------------------------------------------------------------------
    
    Mc = conj(Mc);
    c  = conj(c);
    
    % x component of JIn, store contribution on 3 components of Jout
    fJ = fftn(conj( Mc .* JIn(:,:,:,1)) , [LfK, MfK, NfK]);
    % Jout1 = 0 .* fJ; % First component Vout: 0 * Vin_x
    Jout2 = -fK(:,:,:,3) .* fJ; % Second component Vout: -fK_z * Vin_x
    Jout3 = fK(:,:,:,2) .* fJ; % Third component Vout: +fK_y * Vin_x
    
    % y component of JIn, add contribution on 3 components of Jout
    fJ = fftn(conj( Mc .* JIn(:,:,:,2)) , [LfK, MfK, NfK]);    
    Jout1 = fK(:,:,:,3) .* fJ;  % First component Vout: +fK_z * Vin_y
    % Jout2 = Jout2 + 0 .* fJ; % Second component Vout: 0 * Vin_y
    Jout3 = Jout3 - fK(:,:,:,1) .* fJ; % Third component Vout: -fK_x * Vin_y
    
    % z component of JIn, add contribution on 3 components of Jout
    fJ = fftn(conj( Mc .* JIn(:,:,:,3)) , [LfK, MfK, NfK]);
    Jout1 = Jout1 - fK(:,:,:,2) .* fJ;  % First component Vout: -fK_y * Vin_z
    Jout2 = Jout2 + fK(:,:,:,1) .* fJ; % Second component Vout: +fK_x * Vin_z
    % Jout3 = Jout3 + 0 .* fJ; % Third component Vout: 0 * Vin_z

  
    % apply ifft, multiply by material properties
    Jout1 = ifftn(Jout1);
    JOut(:,:,:,1) = c * conj( Jout1(1:L,1:M,1:N));
    Jout2 = ifftn(Jout2);
    JOut(:,:,:,2) = c *  conj( Jout2(1:L,1:M,1:N));
    Jout3 = ifftn(Jout3);
    JOut(:,:,:,3) = c *  conj( Jout3(1:L,1:M,1:N));
    

elseif strcmp(transp_flag,'notransp') % JOut = fK * JIn0
    
    % ---------------------------------------------------------------------
    % apply fft and mv-op for each of the components of JIn
    % ---------------------------------------------------------------------
       
    % x component of JIn, store contribution on 3 components of Jout
    fJ = fftn(JIn(:,:,:,1),[LfK, MfK, NfK]); 
    % Jout1 = 0 .* fJ; % First component Vout: 0 * Vin_x
    Jout2 = -fK(:,:,:,3) .* fJ; % Second component Vout: -fK_z * Vin_x
    Jout3 = fK(:,:,:,2) .* fJ; % Third component Vout: +fK_y * Vin_x
    
    % y component of JIn, add contribution on 3 components of Jout
    fJ = fftn(JIn(:,:,:,2),[LfK, MfK, NfK]); 
    Jout1 = fK(:,:,:,3) .* fJ;  % First component Vout: +fK_z * Vin_y
    % Jout2 = Jout2 + 0 .* fJ; % Second component Vout: 0 * Vin_y
    Jout3 = Jout3 - fK(:,:,:,1) .* fJ; % Third component Vout: -fK_x * Vin_y
    
    % z component of JIn, add contribution on 3 components of Jout
    fJ = fftn(JIn(:,:,:,3),[LfK, MfK, NfK]); 
    Jout1 = Jout1 - fK(:,:,:,2) .* fJ;  % First component Vout: -fK_y * Vin_z
    Jout2 = Jout2 + fK(:,:,:,1) .* fJ; % Second component Vout: +fK_x * Vin_z
    % Jout3 = Jout3 + 0 .* fJ; % Third component Vout: 0 * Vin_z
    
    % apply ifft, multiply by material properties
    Jout1 = ifftn(Jout1);
    JOut(:,:,:,1) = c .* Mc .* Jout1(1:L,1:M,1:N);
    Jout2 = ifftn(Jout2);
    JOut(:,:,:,2) = c .* Mc .* Jout2(1:L,1:M,1:N);
    Jout3 = ifftn(Jout3);
    JOut(:,:,:,3) = c .* Mc .* Jout3(1:L,1:M,1:N);
    

end

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