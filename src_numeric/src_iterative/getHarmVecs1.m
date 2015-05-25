% getHarmVecs1      For use with GCRODR
%
% Determines harmonic Ritz vectors using matrix H computed from a
% GMRES iteration. For this case, the harmonic Ritz values are the
% eigenvalues of H
% 
% INPUT:  M        dimension of upper Hessenburg matrix H
%         K        select and return basis for space spanned by K harmonic 
%                  Ritz vectors corresponding to K harmonic Ritz values 
%                  of smallest magnitude
%         H        M+1-by-M upper Hessenburg matrix computed from GMRES 
% OUTPUT: HARMVECS basis for span of K harmonic Ritz vectors

function harmVecs = getHarmVecs1(m,k,H)

% Build matrix for eigenvalue problem.
harmRitzMat = H(1:m,:)' \ eye(m);
harmRitzMat(1:m,1:m-1) = 0;
harmRitzMat = H(1:m,:) + H(m+1,m)^2 * harmRitzMat;

% Compute k smallest harmonic Ritz pairs.
[harmVecs, harmVals] = eig(harmRitzMat);
harmVals=abs(diag(harmVals));
% Sort by magnitide of eigenvalue
hv = [harmVals harmVecs'];
hv = sortrows(hv);

% Select k smallest
harmVecs = hv(1:k,2:m+1)';