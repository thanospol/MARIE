function [I_DE] = DIRECT_WS_VA_RWG(r1,r2,r3,r4,r5,Np_1D)
%% Main body of the DIRECT EVALUATION method for the
% evaluation of the vertex adjacent 4-D weakly singular integrals over planar
% triangular elements.

%  Licensing: This code is distributed under the GNU LGPL license. 

%  Modified:  19 October 2011

%  Author:    Athanasios Polimeridis
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Triangle P                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rp1 = r1;
rp2 = r2;
rp3 = r3;
%%%%%%%%%%%%%%%%%%  Triangle Area   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ap = (1/2)*norm(cross(rp2-rp1,rp3-rp1));
%%%%%%%%%%%%%%%%%%  Jacobian   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                
Jp = Ap/sqrt(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Triangle Q                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rq1 = r1;
rq2 = r4;
rq3 = r5;
%%%%%%%%%%%%%%%%%%  Triangle Area   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Aq = (1/2)*norm(cross(rq2-rq1,rq3-rq1));
%%%%%%%%%%%%%%%%%%  Jacobian   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
Jq = Aq/sqrt(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_vo  = (r2-r1)/2;
alpha_v1  = (2*r3-r1-r2)/(2*sqrt(3));
alpha_v2  = (r4-r1)/2;
alpha_v3  = (2*r5-r1-r4)/(2*sqrt(3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[coef_const] = 1.0;
%-------------------------------------------------------------------------%
%                  3-D  Gauss-Legendre Quadrature Rule                    %
%-------------------------------------------------------------------------%
N_theta_p   = Np_1D;
N_theta_q   = Np_1D;
N_psi       = Np_1D;

[w_theta_p,z_theta_p] = Gauss_1D(N_theta_p);
[w_theta_q,z_theta_q] = Gauss_1D(N_theta_q);
[w_psi,z_psi]         = Gauss_1D(N_psi);
%---------------------------------------------
    %%%%%%%%%%%%%%%%%%%%% Gauss quadrature %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    I_theta_p_const = 0;
    
        for n_theta_p = 1:N_theta_p
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%% Int_1, theta_p_A=< theta_p 0 <= theta_p_A  %%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        theta_p_A = 0;
        theta_p_B = pi/3;
        %
        THETA_p = ((theta_p_B-theta_p_A)/2)*z_theta_p(n_theta_p)+(theta_p_B+theta_p_A)/2;
        %
        J_theta_p = (theta_p_B-theta_p_A)/2;
        %
        Lp = (2*sqrt(3))/(sin(THETA_p)+sqrt(3)*cos(THETA_p));
        %
        b_v0 = alpha_vo*cos(THETA_p) + alpha_v1*sin(THETA_p);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%% gauss quadrature %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    I_theta_q_const = 0;
    
        for n_theta_q = 1:N_theta_q
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%% Int_1, theta_q_A=< theta_q 0 <= theta_q_A  %%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        theta_q_A = 0;
        theta_q_B = pi/3;
        %
        THETA_q = ((theta_q_B-theta_q_A)/2)*z_theta_q(n_theta_q)+(theta_q_B+theta_q_A)/2;
        %
        J_theta_q = (theta_q_B-theta_q_A)/2;
        %
        Lq = (2*sqrt(3))/(sin(THETA_q)+sqrt(3)*cos(THETA_q));
        %
        b_v1 = alpha_v2*cos(THETA_q) + alpha_v3*sin(THETA_q);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        a_v = norm(b_v0)^2;
        b_v = -2*dot(b_v0,b_v1);
        c_v = norm(b_v1)^2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        I_psi_const = 0;
        
        for n_psi = 1:N_psi
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%% psi_A =< PSI <= psi_B %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             psi_1_A = 0;
             psi_1_B = atan(Lq/Lp);
             %
             PSI_1 = ((psi_1_B-psi_1_A)/2)*z_psi(n_psi)+(psi_1_B+psi_1_A)/2;
             %
             J_psi_1 = (psi_1_B-psi_1_A)/2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             psi_2_A = atan(Lq/Lp);
             psi_2_B = pi/2;
             %
             PSI_2 = ((psi_2_B-psi_2_A)/2)*z_psi(n_psi)+(psi_2_B+psi_2_A)/2;
             %
             J_psi_2 = (psi_2_B-psi_2_A)/2;
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %
             B_1 =sqrt(a_v*cos(PSI_1)^2+b_v*cos(PSI_1)*sin(PSI_1)+c_v*sin(PSI_1)^2);
             %
             %
             B_2 =sqrt(a_v*cos(PSI_2)^2+b_v*cos(PSI_2)*sin(PSI_2)+c_v*sin(PSI_2)^2);
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %
             Lp = (2*sqrt(3))/(sin(THETA_p)+sqrt(3)*cos(THETA_p));
             L1 = Lp/cos(PSI_1);

             Lq = (2*sqrt(3))/(sin(THETA_q)+sqrt(3)*cos(THETA_q));
             L2 = Lq/sin(PSI_2);
             %
             K_1 = K_functions(B_1,L1);
             
             K_2 = K_functions(B_2,L2);
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %           
             Omega_const_1 = Omega_function_const( PSI_1, B_1, coef_const, K_1);

             Omega_const_2 = Omega_function_const( PSI_2, B_2, coef_const, K_2);
             
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %
             I_psi_const  = I_psi_const  + w_psi(n_psi)*(J_psi_1*Omega_const_1 + J_psi_2*Omega_const_2);
             
        end  %for n_psi = 1:N_psi
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    I_theta_q_const  = I_theta_q_const +w_theta_q(n_theta_q)*I_psi_const;
    
        end  %for n_theta_q = 1:N_theta_q
    %
    I_theta_q_const = J_theta_q*I_theta_q_const;
    %
    I_theta_p_const  = I_theta_p_const +w_theta_p(n_theta_p)*I_theta_q_const;
    
        end  %for n_theta_p = 1:N_theta_p
    %
    I_theta_p_const = J_theta_p*I_theta_p_const;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final Output 
I_DE  = Jp * Jq * I_theta_p_const;
%%
function X = Omega_function_const( Psi, GAMMA, coef, K)

X = K(2) * (sin(Psi) * cos(Psi) * coef(1)) / GAMMA;
%
function K = K_functions(GAMMA,L)

K = zeros(4);

global ko
j = sqrt(-1);

a = j*ko*GAMMA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        K functions                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
K(1) = (1/a^2)*(1-exp(-a*L)-a*L*exp(-a*L));
K(2) = (1/a^3)*(2-2*exp(-a*L)-2*a*L*exp(-a*L)-a^2*L^2*exp(-a*L));
K(3) = -(L^3*exp(-a*L))/a+(3/a)*K(2);
K(4) = -(L^4*exp(-a*L))/a+(4/a)*K(3);