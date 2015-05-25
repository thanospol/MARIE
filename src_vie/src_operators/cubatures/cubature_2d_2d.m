function [I_mn] = cubature_2d_2d(Np_1D,r_m,r_n,R_faces,ko)
%%
N_P = Np_1D;
N_Q = Np_1D;
%
[Np,wp,up,vp] = Gauss_2D(N_P);
[Nq,wq,uq,vq] = Gauss_2D(N_Q);
% N_total = Np*Nq;
I_mn = zeros(6,6);
%
for i_obs = 1:6
    rp1 = R_faces(:,1,i_obs);
    rp2 = R_faces(:,2,i_obs);
    % rp3 = R_faces(:,3,i_obs);
    rp4 = R_faces(:,4,i_obs);
    for i_src = 1:6       
        rq1 = R_faces(:,1,i_src);
        rq2 = R_faces(:,2,i_src);
        % rq3 = R_faces(:,3,i_src);
        rq4 = R_faces(:,4,i_src);
        %                    Opservation face P                            
        a_p = rp2-rp1;
        b_p = rp4-rp1;
        Jp = norm(cross(a_p,b_p));
        %                         Source face P                            
        a_q = rq2-rq1;
        b_q = rq4-rq1;
        Jq = norm(cross(a_q,b_q));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                 4D Numerical Cubature                           %
        %                    Fully vectorized                             % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Weights
        W = wp*wq';
        % Observation points
        Rp_1 = repmat(rp1,1,Np);
        Rm   = repmat(r_m,1,Np);
        rrp = Rp_1 + a_p * up' + b_p * vp' + Rm;
        % Source points
        Rq_1 = repmat(rq1,1,Nq);
        Rn   = repmat(r_n,1,Nq);
        rrq = Rq_1 + a_q * uq' + b_q * vq' + Rn;
        % Distances between observation and source points
        RRR = pdist2(rrp',rrq','euclidean');
        % Singular kernel
        G = exp(-1i*ko*RRR)./ RRR;
        % Integrand
        Integrand =  G;
        % Integral
        Inum = sum(sum(W .* Integrand));
        I_mn(i_obs,i_src) = (Jp*Jq) * Inum;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%