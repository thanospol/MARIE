function [J,E,B,S,Gsar,Pabs] = VWIE_Fields(RHBM,WCOIL,Jc,Jb,freq,flags)
%%    Compute fields due to the solution to the WIE+VIE problem
% _________________________________________________________________________
%
%   Applies operators to generate the fields and other figures of merit
%
% _________________________________________________________________________
%
%% INPUT
%       RHBM  structure
%           r - mapping of the internal edge number to dof number
%           epsilon_r - voxel dielectric
%           sigma_e  - voxel conductivity
%           rho  - voxel density
%           ...
%       WCOIL structure
%           Pcoil - positive end of segment 
%           Ncoil - negative end of segment
%           Dwire - diameter of wire
%           Rhowire - resistivity of material
%           port - port definition
%       Jc - currents in the coil
%       Jb - volumetric currents in the body
%       freq - frequency
%       flags -  for the different computations
%
%
%% OUTPUT
%       J        Solution volumetric currents (LxMxNx3)
%       E        Solution electric field (LxMxNx3)
%       H        Solution magnetic field (LxMxNx3)
%       S        Local SAR (LxMxN)
%       Gsar     Global SAR (number)
%       Pabs     Number
%
%
% -------------------------------------------------------------------------
%
%   J. Fernandez Villena -- jvillena@mit.edu
%   A.G. Polimeridis -- thanos_p@mit.edu
%   Computational Prototyping Group, RLE at MIT
%
% _________________________________________________________________________


if(nargin < 6 )
   flags = squeeze(ones(1,9));
end

tini = tic;

% -------------------------------------------------------------------------
%                 define EM vars and constants
% -------------------------------------------------------------------------

% generate EM constants
EMconstants;

% properties
e_r = RHBM.epsilon_r - 1j*RHBM.sigma_e/(eo*omega);

% dimensions
dx = RHBM.r(2,1,1,1) - RHBM.r(1,1,1,1);
Gram = dx^3;

[L,M,N] = size(e_r);
ND = L*M*N; % number of variables in the system

% get the positions of the non-air voxels in the 3D grid
idxS = find(abs(e_r(:)-1) > 1e-12); % these are the indexes of the non-air voxel positions
idxS3 = [idxS; ND+idxS; 2*ND+idxS]; % the vector of non-air positions for 3 Cartesian components

% x,y,z coordinates of the domain
xd = RHBM.r(:,:,:,1);
yd = RHBM.r(:,:,:,2);
zd = RHBM.r(:,:,:,3);

% Coordinates of all the voxels in the domain
Dcoord = [xd(:), yd(:), zd(:)];

% Coordinates of body in the domain
Scoord = [xd(idxS), yd(idxS), zd(idxS)];

% get the positions of the air voxels in the 3D grid
idxA = setdiff(1:ND,idxS).'; % these are the indexes of the air voxel positions
% Coordinates of air in the domain
Acoord = [xd(idxA), yd(idxA), zd(idxA)];
NA = length(idxA);

clear xd; clear yd; clear zd;

% -------------------------------------------------------------------------
%      prepare the coil current in the required format 
% -------------------------------------------------------------------------

Ccoil = (WCOIL.Pcoil + WCOIL.Ncoil)/2; % center of coil segments
Dcoil = WCOIL.Pcoil - WCOIL.Ncoil; % length of each segment

Jcc = Dcoil; % scale the current by the length of each segment
for ii = 1:3
    Jcc(:,ii) = Jcc(:,ii).*Jc; % for each component
end


% -------------------------------------------------------------------------
%                 Solve the system if needed
% -------------------------------------------------------------------------

t1 = tic;
Einc = []; Hinc = [];
if isempty(Jb)
    
    % -------------------------------------------------------------------------
    %  Generate circulants if needed
   
    if (isempty(RHBM.fN))
        % we need to compute a new circulant
        [RHBM.fN] = getOPERATORS(RHBM.r,freq,'N');
        RHBM.freqfN = freq;
    else
        if (RHBM.freqfN ~= freq)
            [RHBM.fN] = getOPERATORS(RHBM.r,freq,'N');
            RHBM.freqfN = freq;
        end
    end
    
    
    % -------------------------------------------------------------------------
    %  Compute incident fields due to Jc
    
    % [Einc] = E_field_DGF(Jcc,Ccoil,Scoord,ko);
    % [Hinc] = H_field_DGF(Jcc,Ccoil,Scoord,ko);
    [Ebody,Hbody] = eval_DGF(Jcc,Ccoil,Scoord,ko);
    
    % in the whole domain, now reshape them to LMN3 size
    Einc = zeros(L,M,N,3);
    Hinc = zeros(L,M,N,3);
    Einc(idxS3) = Ebody;
    Hinc(idxS3) = Hbody;
    clear Ebody; clear Hbody;
    
    [J,~,~,~,~,~] = JVIE_Solver(Einc,e_r,RHBM.r,freq,RHBM.fN,1e-4);
    
else
    
    J = zeros(L,M,N,3);
    J(idxS3) = Jb;
    
end

fprintf(1,'\n Volumetric currents generated:  \tElapsed time  = %.2f [sec]' , toc(t1));



% -------------------------------------------------------------------------
%                   Compute SAR and GSAR
% -------------------------------------------------------------------------

if flags(4)
    
    t1 = tic;
    
    % get the values of sigma and rho that are reasonable
    denom = 2*RHBM.sigma_e.*RHBM.rho;
    idxsar = find(abs(denom(:))>1e-10);
    idxzero = setdiff(1:ND,idxsar);
    S = J.*conj(J);
    S = sum(S,4);
    S(idxzero) = 0;
    S(idxsar) = S(idxsar)./denom(idxsar);
    Gsar = sum(S(:))*Gram;
    
    % % S= J.*conj(J);
    % % S = sum(S,4);
    % % idxsar = find(abs(RHBM.sigma_e(:)>1e-10));
    % % S(idxsar) = S(idxsar)./(RHBM.sigma_e(idxsar));
    % % idxsar = find(abs(RHBM.rho(:)>1e-12));
    % % S(idxsar) = S(idxsar)./(2*RHBM.rho(idxsar));
    % % Gsar = sum(S(:))*Gram;
    
    fprintf(1,'\n SAR computed:                   \tElapsed time  = %.2f [sec]' , toc(t1));
    
else
    
    S = [];
    Gsar = [];
    
end


% -------------------------------------------------------------------------
%                   Compute Pabs
% -------------------------------------------------------------------------

Pabs = []; Eb = [];
if flags(8)
    
    t1 = tic;

    % -------------------------------------------------------------------------
    %                   Compute E field inside body
    % -------------------------------------------------------------------------

    Eb = J; % copy the currents into the field vector
    Mc = 1j*omega*eo*(e_r-1); % mutiplier to transform currents to fields
    Mc3 = [Mc(:); Mc(:); Mc(:)]; % 3D, for vectors
    Eb(idxS3) = Eb(idxS3)./Mc3(idxS3); % transform the non-air positions
    % note that it will be zero outside
    
    % ---------------------------------------------------------------------
    %                   Compute Pabs in body
    % ---------------------------------------------------------------------
    
    Pabs = Eb.*conj(J);
    Pabs = sum(Pabs,4);
    Pabs = 0.5*real(sum(Pabs(:)))*Gram;    
    
    fprintf(1,'\n Absorbed Power computed:        \tElapsed time  = %.2f [sec]' , toc(t1));
    
end



% -------------------------------------------------------------------------
%                   Compute E field
% -------------------------------------------------------------------------

E = [];
if flags(5)
    
    t1 = tic;
    
    if isempty(Eb)
        
        % -------------------------------------------------------------------------
        %                   Compute E field inside body
        % -------------------------------------------------------------------------
        
        Eb = J; % copy the currents into the field vector
        Mc = 1j*omega*eo*(e_r-1); % mutiplier to transform currents to fields
        Mc3 = [Mc(:); Mc(:); Mc(:)]; % 3D, for vectors
        Eb(idxS3) = Eb(idxS3)./Mc3(idxS3); % transform the non-air positions
        % note that it will be zero outside
        
    end
    
    if flags(9) % need E field in the whole domain
        
        if isempty(Einc) % we need the incident field and circulant
            
            % ---------------------------------------------------------------------
            %                  Generate circulants if needed
            if (isempty(RHBM.fN))
                % we need to compute a new circulant
                [RHBM.fN] = getOPERATORS(RHBM.r,freq,'N');
                RHBM.freqfN = freq;
            else
                if (RHBM.freqfN ~= freq)
                    [RHBM.fN] = getOPERATORS(RHBM.r,freq,'N');
                    RHBM.freqfN = freq;
                end
            end
            
            
            % we only compute the E fields outside the body... in body is ready
            
            % -----------------------------------------------------------------
            %         Compute incident fields due to Jc
            
            % [Einc] = E_field_DGF(Jcc,Ccoil,Scoord,ko);
            % [Hinc] = H_field_DGF(Jcc,Ccoil,Scoord,ko);
            if flags(6)
                % E and H field, point matching in whole domain
                [Einc,Hinc] = eval_DGF(Jcc,Ccoil,Dcoord,ko);
                % in the whole domain, now reshape them to LMN3 size
                Einc = reshape(Einc,L,M,N,3);
                Hinc = reshape(Hinc,L,M,N,3);
                
            else
                % only E field, and we only need in the Air
                [EincAir] = E_field_DGF(Jcc,Ccoil,Acoord,ko);
                Einc = zeros(L,M,N,3);
                EincComp = zeros(L,M,N); % component
                EincAir = reshape(EincAir,NA,3);
                for ii = 1:3
                    EincComp(idxA) = EincAir(:,ii);
                    Einc(:,:,:,ii) = EincComp;
                end
                clear EincAir; clear EincComp;
                
            end
            
            % Note that in general we only have valid Einc in air positions
            
        end
        
        % ---------------------------------------------------------------------
        %                   Compute E field
        % ---------------------------------------------------------------------
        
        [E] = E_field_Nop(J,RHBM.fN,Gram,freq,Einc);
        
        % overwrite E in non-air voxels with the ones we already have
        E(idxS3) = Eb(idxS3);
        
    else % we only need E field in the body
        
        E = Eb; % we already had them
        
    end
    
    fprintf(1,'\n E fields computed:              \tElapsed time  = %.2f [sec]' , toc(t1));
    
end

clear Eb;

        
% -------------------------------------------------------------------------
%                   Compute H field
% -------------------------------------------------------------------------

B = [];
if flags(6)
    
    t1 = tic;
    
    if isempty(Hinc)
        
        % -----------------------------------------------------------------
        %         Compute incident fields due to Jc
        
        if flags(9) % we need the fields in the whole domain
            
            % Point matching for the H field in Dcoord
            [Hinc] = H_field_DGF(Jcc,Ccoil,Dcoord,ko);
            % in the whole domain, now reshape them to LMN3 size
            Hinc = reshape(Hinc,L,M,N,3);
            
        else % we only need fields in the body
            
            % Point matching for the H field in Scoord
            [Hbody] = H_field_DGF(Jcc,Ccoil,Scoord,ko);
            Hinc = zeros(L,M,N,3);
            Hinc(idxS3) = Hbody;
            clear Hbody;
            
        end
        
    end
        
    
    % ---------------------------------------------------------------------
    %        Generate circulants if needed
    
    if (isempty(RHBM.fK))
        % we need to compute a new circulant
        [RHBM.fK] = getOPERATORS(RHBM.r,freq,'K');
        RHBM.freqfK = freq;
    else
        if (RHBM.freqfK ~= freq)
            [RHBM.fK] = getOPERATORS(RHBM.r,freq,'K');
            RHBM.freqfK = freq;
        end
    end
    

    % ---------------------------------------------------------------------
    %         Compute Magnetic Field
    % ---------------------------------------------------------------------
        
    [Ball] = H_field_Kop(J,RHBM.fK,Gram,Hinc);
     
    if flags(9) % fields in the whole domain
        B = 4*pi*1e-7*Ball;
    else % we only need fields in the body
        % Ball has the scattered field in the air positions: need to zero it
        B = zeros(L,M,N,3);
        B(idxS3) = 4*pi*1e-7*Ball(idxS3);   
    end
    
    fprintf(1,'\n B fields computed:              \tElapsed time  = %.2f [sec]' , toc(t1));

    
end


% -------------------------------------------------------------------------
%                 And it is done
% -------------------------------------------------------------------------

toverall = toc(tini);
fprintf(1,'\n --------------------------------------------------------------');
fprintf(1,'\n Overall:                        \tElapsed time  = %.2f [sec]', toverall);
