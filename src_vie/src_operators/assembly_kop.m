function [IK_mn] = assembly_kop(r,ko,dx)
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Grid Dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[L, M, N, ~] = size(r);
dy =dx; dz = dx;

% Get the 6 faces of  the cube
[R_faces] = cube_faces(dx,dy,dz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Quadrature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Np_1D_near   = 20;
Np_1D_medium = 10;
Np_1D_far    = 4;
% Np_GL = [Np_DEMCEM ; Np_1D_near ; Np_1D_medium ; Np_1D_far];

% Set region of medium distances cells for higher order quadrature
n_medium = 5; 

% Allocate memory for main matrices
IK_mn = zeros(L,M,N,3);
% Reference cell
r_n = [0.0 , 0.0 , 0.0]';
% Reference distance vector
d = [dx,dy,dz];
% modification for "parfor"
RR = R_faces;
kko = ko;

tic_Assembly = tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            ASSEMBLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Far Distance Cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_P = Np_1D_far;

[Np,wp,up,vp] = gauss_2d(N_P);

parfor mx = 1:L
    for my = 1:M
        for mz = 1:N
            m = [mx,my,mz];
            r_m = ( (m-1) .* d )';

            IK_mn(mx,my,mz,:)  = mexcubature_kop(Np,wp,up,vp,r_m,r_n,RR,kko);
                     
        end
    end
end

Time_far = toc;
fprintf(' Time_far         = %d \n',int64(Time_far));
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Medium Distance Cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_P = Np_1D_medium;

[Np,wp,up,vp] = gauss_2d(N_P);

parfor mx = 1:n_medium
    for my = 1:n_medium
        for mz = 1:n_medium
            
            m = [mx,my,mz];
            r_m = ( (m-1) .* d )';

            IK_mn(mx,my,mz,:)  = mexcubature_kop(Np,wp,up,vp,r_m,r_n,RR,kko);

                     
        end
    end
end

Time_medium = toc;
fprintf(' Time_medium      = %d \n',int64(Time_medium));
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Nearby Cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_P = Np_1D_near;

[Np,wp,up,vp] = gauss_2d(N_P);

parfor mx = 1:2
    for my = 1:2
        for mz = 1:2
            
            m = [mx,my,mz];
            r_m = ( (m-1) .* d )';
                       
            IK_mn(mx,my,mz,:)  = cubature_kop_near_sym(Np,wp,up,vp,Np,wp,up,vp,r_m,r_n,m,RR,kko);
 
        end
    end
end

Time_near = toc;
fprintf(' Time_near        = %d \n',int64(Time_near));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     FINAL OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coef = dx^4 / (1i*ko)^2 / (4*pi);
%
IK_mn = coef * IK_mn;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             END ASSEMBLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Time_Assembly = toc(tic_Assembly);
fprintf(' ------------------------  \n')
fprintf(' Time_Assembly    = %dm%ds  \n' ,floor(Time_Assembly/60),int64(mod(Time_Assembly,60)))