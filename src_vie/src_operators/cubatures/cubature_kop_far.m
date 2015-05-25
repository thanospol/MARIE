function [IK_mn] = cubature_kop_far(Np,wp,up,vp,Nq,wq,uq,vq,r_m,r_n,RR,kko)
%%
% based on E. H. Bleszynski et al. IEEE-TAP 2013 
R_faces = RR;
ko = kko;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
index_obs    = [1  1  2  2  3  3  4  4  5  5  6  6]';
index_src    = [1  2  1  2  3  4  3  4  5  6  5  6]';
sign_obs_src = [1 -1 -1  1  1 -1 -1  1  1 -1 -1  1]';
%
% N_unique = length(index_obs);
N_unique = 12;
%
I_mn_x = zeros(N_unique,1);
I_mn_y = zeros(N_unique,1);
I_mn_z = zeros(N_unique,1);
%
Rm   = repmat(r_m,1,Np);
Rn   = repmat(r_n,1,Nq);
%
Rmn_x = zeros(Np,Np);
Rmn_y = zeros(Np,Np);
Rmn_z = zeros(Np,Np);
%%
for ind_unique = 1:N_unique  
    % Observation Panel
    i_obs = index_obs(ind_unique);
    %
    rp1 = R_faces(:,1,i_obs);
    rp2 = R_faces(:,2,i_obs);
    % rp3 = R_faces(:,3,i_obs);
    rp4 = R_faces(:,4,i_obs);
    %              
    a_p = rp2-rp1;
    b_p = rp4-rp1;
    % Source Panel
    i_src = index_src(ind_unique);
    %
    rq1 = R_faces(:,1,i_src);
    rq2 = R_faces(:,2,i_src);
    % rq3 = R_faces(:,3,i_src);
    rq4 = R_faces(:,4,i_src);
    %         
    a_q = rq2-rq1;
    b_q = rq4-rq1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                 4D Numerical Cubature                           %
    %                    Fully vectorized                             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Weights
    W = wp*wq';
    % Observation points
    Rp_1 = repmat(rp1,1,Np);
    rrp = Rp_1 + a_p * up' + b_p * vp' + Rm;
    % Source points
    Rq_1 = repmat(rq1,1,Nq);
    rrq = Rq_1 + a_q * uq' + b_q * vq' + Rn;
    % Distances between observation and source points instead of Rpq = pdist2(rrp',rrq','euclidean');
    for ii=1:Np
        Rmn_x(ii,:) = rrp(1,ii)-rrq(1,:);
        Rmn_y(ii,:) = rrp(2,ii)-rrq(2,:);
        Rmn_z(ii,:) = rrp(3,ii)-rrq(3,:);
    end
    %
    Rpq = sqrt(Rmn_x.^2 + Rmn_y.^2 + Rmn_z.^2);
    % Singular kernel
    G = (1-exp(-1i*ko*Rpq) .* (1 + 1i*ko*Rpq) )./ Rpq.^3;
    % Integral
    I_mn_x(ind_unique) =  sum(sum(W .* G .* Rmn_x));
    I_mn_y(ind_unique) =  sum(sum(W .* G .* Rmn_y));
    I_mn_z(ind_unique) =  sum(sum(W .* G .* Rmn_z));
end
%
IK_mn_x = sum(I_mn_x.*sign_obs_src);
IK_mn_y = sum(I_mn_y.*sign_obs_src);
IK_mn_z = sum(I_mn_z.*sign_obs_src);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          FINAL OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IK_mn = [IK_mn_x IK_mn_y IK_mn_z]'; % dx^2 = area of panel = Jacobian