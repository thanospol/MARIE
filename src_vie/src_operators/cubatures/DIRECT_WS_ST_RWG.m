function [I_DE] = DIRECT_WS_ST_RWG(r1,r2,r3,Np_1D)
%% Main body of the DIRECT EVALUATION method for the
% evaluation of the coincident 4-D weakly singular integrals over planar
% triangular elements.

%  Licensing: This code is distributed under the GNU LGPL license. 

%  Modified:  19 October 2011

%  Author:    Athanasios Polimeridis

% References

% A. G. Polimeridis and T. V. Yioultsis, “On the direct evaluation of weakly singular
% integrals in Galerkin mixed potential integral equation formulations,” IEEE Trans.
% Antennas Propag., vol. 56, no. 9, pp. 3011-3019, Sep. 2008.

% A. G. Polimeridis and J. R. Mosig, “Complete semi-analytical treatment of weakly
% singular integrals on planar triangles via the direct evaluation method,” Int. J.
% Numerical Methods Eng., vol. 83, pp. 1625-1650, 2010.

% A. G. Polimeridis, J. M. Tamayo, J. M. Rius and J. R. Mosig, “Fast and accurate
% computation of hyper-singular integrals in Galerkin surface integral equation
% formulations via the direct evaluation method,” IEEE Trans.
% Antennas Propag., vol. 59, no. 6, pp. 2329-2340, Jun. 2011.

% A. G. Polimeridis and J. R. Mosig, “On the direct evaluation of surface integral
% equation impedance matrix elements involving point singularities,” IEEE Antennas
% Wireless Propag. Lett., vol. 10, pp. 599-602, 2011.

% INPUT DATA
% r1,r2,r3 = point vectors of the triangular element's vertices
% Outer triangle P:(rp1,rp2,rp3)=(r1,r2,r3)
% Inner triangle Q:(rq1,rq2,rq3)=(r1,r2,r3)
% Np_1D = order of the Gauss-Legendre quadrature rule

% OUTPUT DATA
%             I_DE(1)  = I_f1_f1
%             I_DE(2)  = I_f1_f2
%             I_DE(3)  = I_f1_f3
%             I_DE(4)  = I_f2_f1
%             I_DE(5)  = I_f2_f2
%             I_DE(6)  = I_f2_f3
%             I_DE(7)  = I_f3_f1
%             I_DE(8)  = I_f3_f2
%             I_DE(9)  = I_f3_f3
%
%%
%%%%%%%%%%%%%%%%%%  Triangle Area   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ap = (1/2)*norm(cross(r2-r1,r3-r1));
%%%%%%%%%%%%%%%%%%  Jacobian   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                
Jp = Ap/sqrt(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Asing0 = norm(r1-r2)^2/4;
Asing1 = dot(r1-r2,r1+r2-2*r3)/(2*sqrt(3));
Asing2 = norm(r1+r2-2*r3)^2/12; 
%
Asing3 = (1/4)        *Asing0   -(sqrt(3)/4)  *Asing1      +(3/4)        *Asing2;
Asing4 = (sqrt(3)/2)  *Asing0   -(1/2)        *Asing1      -(sqrt(3)/2)  *Asing2;
Asing5 = (3/4)        *Asing0   +(sqrt(3)/4)  *Asing1      +(1/4)        *Asing2;
%
Asing6 = (1/4)        *Asing0   +(sqrt(3)/4)  *Asing1       +(3/4)       *Asing2;
Asing7 = -(sqrt(3)/2) *Asing0   -(1/2)        *Asing1       +(sqrt(3)/2) *Asing2;
Asing8 = (3/4)        *Asing0   -(sqrt(3)/4)  *Asing1       +(1/4)       *Asing2;
%
Asing = [Asing0 Asing3 Asing6
         Asing1 Asing4 Asing7
         Asing2 Asing5 Asing8];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  1D  Gauss-Legendre Quadrature Rule                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[w,z] = Gauss_1D(Np_1D);

Isub = zeros(3,3,3);
%---------------------------------------------
for m = 1:3
    %%%%%%%%% Define the coefs of the appropriate subtriangle
    acc = Asing(1,m);
    acs = Asing(2,m);
    ass = Asing(3,m);
    %%%%%%%%%%%%%%%%%%%%% Gauss quadrature %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    Int1_1_1=0;Int2_1_1=0;Int3_1_1=0;
    for kk = 1:Np_1D
        %%%%%%%%%%%%%%% Int1,  0 =< PSI <= pi/3 %%%%%%%%%%%%%%%%%%%%%%%%%%%
        PSI_1a = 0;
        PSI_1b = pi/3;
        PSI_1k = ((PSI_1b-PSI_1a)*z(kk)+(PSI_1b+PSI_1a))/2;
        a_PSI_1k = sqrt(acc*cos(PSI_1k)^2-acs*cos(PSI_1k)*sin(PSI_1k)+ass*sin(PSI_1k)^2);
        % Int1_1_1 
        PHI_a_1_1 = Phi_1_1(1,PSI_1k,a_PSI_1k);
        PHI_e_1_1 = Phi_1_1(5,PSI_1k,a_PSI_1k);
        PHI_f_1_1 = Phi_1_1(6,PSI_1k,a_PSI_1k);
        
        F1_1_1 = PHI_a_1_1+PHI_e_1_1+PHI_f_1_1;       
        Int1_1_1  = Int1_1_1+w(kk)*F1_1_1;
        %%%%%%%%%%%%%%% Int2,  pi/3 =< PSI <= 2pi/3 %%%%%%%%%%%%%%%%%%%%%%%
        PSI_2a = pi/3;
        PSI_2b = 2*pi/3;
        PSI_2k = ((PSI_2b-PSI_2a)*z(kk)+(PSI_2b+PSI_2a))/2;
        a_PSI_2k = sqrt(acc*cos(PSI_2k)^2-acs*cos(PSI_2k)*sin(PSI_2k)+ass*sin(PSI_2k)^2);
        % Int2_1_1 
        PHI_b_1_1 = Phi_1_1(2,PSI_2k,a_PSI_2k);
        PHI_g_1_1 = Phi_1_1(7,PSI_2k,a_PSI_2k);
                
        F2_1_1 = PHI_b_1_1+PHI_g_1_1;       
        Int2_1_1  = Int2_1_1+w(kk)*F2_1_1;
        %%%%%%%%%%%%%%% Int3,  2pi/3 =< PSI <= pi %%%%%%%%%%%%%%%%%%%%%%%%%
        PSI_3a = 2*pi/3;
        PSI_3b = pi;
        PSI_3k = ((PSI_3b-PSI_3a)*z(kk)+(PSI_3b+PSI_3a))/2;
        a_PSI_3k = sqrt(acc*cos(PSI_3k)^2-acs*cos(PSI_3k)*sin(PSI_3k)+ass*sin(PSI_3k)^2);
        % Int3_1_1 
        PHI_c_1_1 = Phi_1_1(3,PSI_3k,a_PSI_3k);
        PHI_d_1_1 = Phi_1_1(4,PSI_3k,a_PSI_3k);
        PHI_h_1_1 = Phi_1_1(8,PSI_3k,a_PSI_3k);
                
        F3_1_1 = PHI_c_1_1+PHI_d_1_1+PHI_h_1_1;       
        Int3_1_1  = Int3_1_1+w(kk)*F3_1_1;
    end % for kk=1:Npp
    % Isub_1_1
     Int1_1_1 = ((PSI_1b-PSI_1a)/2)*Int1_1_1;
     Int2_1_1 = ((PSI_2b-PSI_2a)/2)*Int2_1_1;
     Int3_1_1 = ((PSI_3b-PSI_3a)/2)*Int3_1_1;
          
     Isub(1,1,m) = Int1_1_1+Int2_1_1+Int3_1_1;
     
end % for m = 1:3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_DE = Jp^2 * (Isub(1,1,1)+Isub(1,1,2)+Isub(1,1,3));

function [PHI] = Phi_1_1(argument,psi,Apsi)
%% Phi_1_1 function

%  Licensing: This code is distributed under the GNU LGPL license. 

%  Modified:  20 September 2011

%  Author:    Athanasios Polimeridis

% References

% A. G. Polimeridis and T. V. Yioultsis, “On the direct evaluation of weakly singular
% integrals in Galerkin mixed potential integral equation formulations,” IEEE Trans.
% Antennas Propag., vol. 56, no. 9, pp. 3011-3019, Sep. 2008.

% A. G. Polimeridis and J. R. Mosig, “Complete semi-analytical treatment of weakly
% singular integrals on planar triangles via the direct evaluation method,” Int. J.
% Numerical Methods Eng., vol. 83, pp. 1625-1650, 2010.

% A. G. Polimeridis, J. M. Tamayo, J. M. Rius and J. R. Mosig, “Fast and accurate
% computation of hyper-singular integrals in Galerkin surface integral equation
% formulations via the direct evaluation method,” IEEE Trans.
% Antennas Propag., vol. 59, no. 6, pp. 2329-2340, Jun. 2011.

% A. G. Polimeridis and J. R. Mosig, “On the direct evaluation of surface integral
% equation impedance matrix elements involving point singularities,” IEEE Antennas
% Wireless Propag. Lett., vol. 10, pp. 599-602, 2011.

% INPUT DATA
% argument  = 1-8 -> a-h
% psi       = Variable Psi
% Apsi      = alpha(Psi)

% OUTPUT DATA
% Phi_1_1   
%
global ko
j = sqrt(-1);

A = j*ko*Apsi;
B = cos(psi);
C = sin(psi)/sqrt(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Int_1_1                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
D = sin(psi)/((j*ko)^2*Apsi^3);

switch argument

case 1
%%%%%%%%%%%%%%%%%%%%%%    PHI_a   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI_a1 = A/(2*B);
PSI_a2 = -1;
PSI_a3 = (B/A)*(1-exp(-A/B));
PSI_a  = PSI_a1+PSI_a2+PSI_a3;
PHI    = D*PSI_a;

case 2
%%%%%%%%%%%%%%%%%%%%%%    PHI_b   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI_b1 = A/(2*C);
PSI_b2 = -1;
PSI_b3 = (C/A)*(1-exp(-A/C));
PSI_b  = PSI_b1+PSI_b2+PSI_b3;
PHI    = D*PSI_b;

case 3
%%%%%%%%%%%%%%%%%%%%%%    PHI_c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta    = tan(pi-psi)/sqrt(3);
gamma   = (1-beta)/(1+beta);

PSI_c1 = (A/C)*(1/2-(gamma-gamma^2/2));
PSI_c2 = -(1-gamma);
PSI_c3 = (C/A)*(1-exp(-A*(1-gamma)/C));
PSI_c  = PSI_c1+PSI_c2+PSI_c3;
PHI    = D*PSI_c;

case 4
%%%%%%%%%%%%%%%%%%%%%%    PHI_d   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta    = tan(pi-psi)/sqrt(3);
gamma   = (1-beta)/(1+beta);

PSI_d1 = -A*(gamma+gamma^2/2)/B;
PSI_d2 = -gamma;
PSI_d3 = (B/A)*(exp(A*(1+gamma)/B)-exp(A/B));
PSI_d  = PSI_d1+PSI_d2+PSI_d3;
PHI    = D*PSI_d;

case 5
%%%%%%%%%%%%%%%%%%%%%%    PHI_e   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon = tan(psi)/sqrt(3);
delta   = -(1-epsilon)/(1+epsilon);

PSI_e1 = A*(delta^2/2-delta)/B;
PSI_e2 = delta;
PSI_e3 = (B/A)*(exp(-A/B)-exp(-A*(1-delta)/B));
PSI_e  = PSI_e1+PSI_e2+PSI_e3;
PHI    = D*PSI_e;

case 6
%%%%%%%%%%%%%%%%%%%%%%    PHI_f   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon = tan(psi)/sqrt(3);
delta   = -(1-epsilon)/(1+epsilon);

PSI_f1 = (A/C)*(delta+delta^2/2+1/2);
PSI_f2 = -(1+delta);
PSI_f3 = (C/A)*(1-exp(-A*(1+delta)/C));
PSI_f  = PSI_f1+PSI_f2+PSI_f3;
PHI    = D*PSI_f;

case 7
%%%%%%%%%%%%%%%%%%%%%%    PHI_g   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI_g1 = A/(2*C);
PSI_g2 = -1;
PSI_g3 = (C/A)*(1-exp(-A/C));
PSI_g  = PSI_g1+PSI_g2+PSI_g3;
PHI    = D*PSI_g;

case 8
%%%%%%%%%%%%%%%%%%%%%%    PHI_h   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI_h1 = -A/(2*B);
PSI_h2 = -1;
PSI_h3 = (B/A)*(exp(A/B)-1);
PSI_h  = PSI_h1+PSI_h2+PSI_h3;
PHI    = D*PSI_h;

end %switch argument
