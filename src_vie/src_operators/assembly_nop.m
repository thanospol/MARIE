function [G_mn] = assembly_nop(r,ko,dx,singular)
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
Np_DEMCEM    = 20;  
Np_1D_near   = 20;
Np_1D_medium = 10;
Np_1D_far    = 4;
% Np_GL = [Np_DEMCEM ; Np_1D_near ; Np_1D_medium ; Np_1D_far];

% Set region of medium distances cells for higher order quadrature
n_medium = 5; 

% Allocate memory for main matrices
J_mn = zeros(L,M,N,6,6);
G_mn = zeros(L,M,N,6);
% Reference cell
r_n = [0.0 , 0.0 , 0.0]';
% Reference distance vector
d = [dx,dy,dz];

tic_Assembly = tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            ASSEMBLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_P = Np_1D_far;

[Np,wp,up,vp] = gauss_2d(N_P);

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Far Distance Cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RR = R_faces;
kko = ko;
ddx = dx;
%
parfor mx = 1:L
    for my = 1:M
        for mz = 1:N
            m = [mx,my,mz];
            r_m = ( (m-1) .* d )';

            J_mn(mx,my,mz,:,:) = mexcubature_nop(Np,wp,up,vp,r_m,r_n,RR,kko,ddx);
                     
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

            J_mn(mx,my,mz,:,:) = mexcubature_nop(Np,wp,up,vp,r_m,r_n,RR,kko,ddx);

                     
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

% Choose singular integration method
[I_ST]        = singular_st(Np_DEMCEM,ko,dx,dy,singular);
[I_EAc,I_EAo] = singular_ea(Np_DEMCEM,ko,dx,dy,dz,singular);
[I_VAc,I_VAo] = singular_va(Np_DEMCEM,ko,dx,dy,dz,singular);

parfor mx = 1:2
    for my = 1:2
        for mz = 1:2
            
            m = [mx,my,mz];
            r_m = ( (m-1) .* d )';
      
%             [I_mn_non_corrected] = cubature_2d_2d(Np_1D_near,r_m,r_n,RR,kko);
            [I_mn_non_corrected] = cubature_nop_sym_near(Np,wp,up,vp,r_m,r_n,RR,kko,ddx);
            
            J_mn(mx,my,mz,:,:) = nearby_correct(m,I_mn_non_corrected,I_ST,I_EAc,I_EAo,I_VAc,I_VAo);
 
        end
    end
end

Time_near = toc;
fprintf(' Time_near        = %d \n',int64(Time_near));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Evaluate G_mn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
tic

[P_vector] = vector_projection();

for kk = 1:6
    for jj = 1:6
        for ii = 1:6
            G_mn(:,:,:,kk) = G_mn(:,:,:,kk) + J_mn(:,:,:,ii,jj)*P_vector(ii,jj,kk);
        end
    end
end

Time_projection = toc;
fprintf(' Time_projection  = %d \n',int64(Time_projection));

G_mn = ( 1.0 / (4*pi) ) * G_mn;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             END ASSEMBLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Time_Assembly = toc(tic_Assembly);
fprintf(' ------------------------  \n')
fprintf(' Time_Assembly    = %dm%ds  \n' ,floor(Time_Assembly/60),int64(mod(Time_Assembly,60)))
