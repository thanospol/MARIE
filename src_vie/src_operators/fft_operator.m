function [fA] = fft_operator(A)
% 
[L,M,N,d] = size(A);
fA = zeros(L,M,N,d);
% 3D-FFT of A
for p=1:d
    fA(:,:,:,p) = fftn(A(:,:,:,p));
end