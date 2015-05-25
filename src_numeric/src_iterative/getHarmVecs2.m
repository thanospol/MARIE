% getHarmVecs2     For use with GCRODR
%
% Determines harmonic Ritz vectors using matrices computed from
% GMRES iteration. 
% 
% INPUT:  M        dimension of upper Hessenburg matrix H
%         K        select and return basis for space spanned by K harmonic 
%                  Ritz vectors corresponding to K harmonic Ritz values 
%                  of smallest magnitude
%         H2       upper Hessenburg matrix computed GCRODR relations
%         V        N-by-M+1 matrix containing Krylov basis computed by GMRES
%         U        basis for recycled subspace
%         C        C = A*U, as per GCRODR relations
% OUTPUT: HARMVECS basis for span of K harmonic Ritz vectors
function harmVecs = getHarmVecs2(m,k,H2,V,U,C)

B = H2' * H2;

% A = | C'*U        0 |
%     | V_{m+1}'*U  I |
A = zeros(m+1,m);
A(1:k,1:k) = C' * U;
A(k+1:m+1,1:k) = V' * U;
A(k+1:m,k+1:m) = eye(m-k);
A = H2' * A;

% Compute k smallest harmonic Ritz pairs.
[harmVecs, harmVals] = eig(A,B);
% Sort by magnitide of eigenvalue
harmVals = abs(diag(harmVals));
hv = [harmVals harmVecs'];
hv = sortrows(hv);

% k smallest harmonic ritz values
% Actually, k largest of (1./harmonic Ritz value)
harmVecs = hv((m-k)+1:m,2:m+1)';