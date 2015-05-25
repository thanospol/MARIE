function [ZP,Jc,Jb,Sb,Eb,Bb,Gsar,Pabs] = MR_Solver(RHBM,COIL,freq,tol,flags,logfile)
%%    Simple driver to completely solve MR EM problem
% _________________________________________________________________________
%
%       Completely solves the MR EM problem with coupled formulation
%
% _________________________________________________________________________
%
%% INPUT
%   RHBM       Body model structure
%   COIL       Coil model structure
%   freq       Frequency vector (in Hz)
%
%
%% OPTIONAL INPUT
%   rho         Density (LxMxN)
%   tol         Relative tolerance for the method (DEFAULT 1e-3)
%   flags       1 or 0 for the selection of the output
%
%
%% OUTPUT
%   ZP          Z parameter matrix (NpxNpxNfreqs)
%   Jc          Coil basis (current) coefficients (NcxNpxNfreqs)
%   Jb          Body basis (current) coefficients (LxMxNx3xNpxNfreqs)
%   Sb          Local SAR (LxMxNxNpxNfreqs)
%   Eb          Solution electric field (LxMxNx3xNpxNfreqs)
%   Bb          Solution magnetic field (LxMxNx3xNpxNfreqs)
%   Gsar        Global SAR (number) (NpxNfreqs)
%   Pabs        Absorbed power (number) (NpxNfreqs)
%
%
% -------------------------------------------------------------------------
%
%   J. Fernandez Villena -- jvillena@mit.edu
%   A.G. Polimeridis -- thanos_p@mit.edu
%   Computational Prototyping Group, RLE at MIT
%
% _________________________________________________________________________



% -------------------------------------------------------------------------
% Initialization of variables
% -------------------------------------------------------------------------

if(nargin < 3 )
   fprintf(1, '\n ERROR: not enough arguments\n');
   return
end
if(nargin < 4 || isempty(tol))
   tol = 1e-3;
else
    if(tol < 1e-15) || (tol >= 1)
        tol = 1e-3;
    end
end
if(nargin < 5 || isempty(flags))
   flags = ones(9,1); % returns all
end

if(nargin < 6 || isempty(logfile))
   fid = 1;
   logfile = [];
else
    fid = fopen(logfile, 'a');
    if (fid == -1)
        fid = 1; % standard output
    end
end

ZP = [];
Jc = [];
Jb = [];
Sb = [];
Eb = [];
Bb = [];
Gsar = [];
Pabs = [];


% -------------------------------------------------------------------------
%   check the existence of different cases
% -------------------------------------------------------------------------

freq = squeeze(freq);
Nfreqs = length(freq); % number of points
[L,M,N,~] = size(RHBM.r);

body = nnz(abs(RHBM.epsilon_r(:) - 1 + 1j*RHBM.sigma_e(:))); % check if there is scatterer

% -------------------------------------------------------------------------
if body % there is scatterer
    
    switch COIL.type

        % -----------------------------------------------------------------
        case 'W' % wire coil model
        
            Nports = length(COIL.port); % get number of ports
            Nvars = size(COIL.Pcoil,1); % number of segments

            % allocate space
            if flags(1)
                ZP = zeros(Nports,Nports,Nfreqs);
            end
            if flags(2)
                Jc = zeros(Nvars,Nports,Nfreqs);
            end
            if flags(3)
                Jb = zeros(L,M,N,3,Nports,Nfreqs);
            end
            if flags(4)
                Sb = zeros(L,M,N,Nports,Nfreqs);
            end
            if flags(5)
                Eb = zeros(L,M,N,3,Nports,Nfreqs);
            end
            if flags(6)
                Bb = zeros(L,M,N,3,Nports,Nfreqs);
            end
            if flags(7)
                Gsar = zeros(Nports,Nfreqs);
            end
            if flags(8)
                Pabs = zeros(Nports,Nfreqs);
            end
            
            for ii = 1:length(freq) % loop on the frequencies

                ff = freq(ii);
                
                tini = tic;
                fprintf(1,'\n\n ---------------------------------------------------------------------');
                fprintf(1,'\n ---------------------------------------------------------------------');
                fprintf(1,'\n Solving the WIE+VIE system at freq %.3f MHz\n ', ff/1e6);
                

                coup = 5; gpu_flag = 1;
                factorized = 1; tolfactor = 1e-5; blockSize = 100;
                % solve the system
                [Zparam,Jcoil,Jbody] = VWIE_Solver(RHBM,COIL,ff,tol,coup,gpu_flag,factorized,tolfactor,blockSize);

                % store solution
                if flags(1)
                    ZP(:,:,ii) = Zparam;
                end
                if flags(2)
                    for jj = 1:Nports
                        Jc(:,jj,ii) = Jcoil(:,jj);
                    end
                end
                if sum(flags(3:8))
                    
                    if (flags(5) && flags(9)) % need electric field in air
                        % -------------------------------------------------
                        %   Generate circulants if needed
                        if (isempty(RHBM.fN))
                            % we need to compute a new circulant
                            [RHBM.fN] = getOPERATORS(RHBM.r,ff,'N');
                            RHBM.freqfN = ff;
                        else
                            if (RHBM.freqfN ~= ff)
                                [RHBM.fN] = getOPERATORS(RHBM.r,ff,'N');
                                RHBM.freqfN = ff;
                            end
                        end
                    end
                    
                    if flags(6)
                        % -------------------------------------------------
                        %                  Generate circulants if needed
                        if (isempty(RHBM.fK))
                            % we need to compute a new circulant
                            [RHBM.fK] = getOPERATORS(RHBM.r,ff,'K');
                            RHBM.freqfK = ff;
                        else
                            if (RHBM.freqfK ~= ff)
                                [RHBM.fK] = getOPERATORS(RHBM.r,ff,'K');
                                RHBM.freqfK = ff;
                            end
                        end
                    end
                    

                    for jj = 1:Nports
                        
                        fprintf(1,'\n\n ----------------');
                        fprintf(1,'\n Port # %d', jj);
                        
                        % compute fields for each channel
                        [J,E,B,S,G,P] = VWIE_Fields(RHBM,COIL,Jcoil(:,jj),Jbody(:,jj),ff,flags);

                        if flags(3)
                            Jb(:,:,:,:,jj,ii) = J;
                        end
                        if flags(4)
                            Sb(:,:,:,jj,ii) = S;
                        end
                        if flags(5)
                            Eb(:,:,:,:,jj,ii) = E;
                        end
                        if flags(6)
                            Bb(:,:,:,:,jj,ii) = B;
                        end
                        if flags(7)
                            Gsar(jj,ii) = G;
                        end
                        if flags(8)
                            Pabs(jj,ii) = P;
                        end

                    end

                end
                
                fprintf(1,'\n\n ---------------------------------------------------------------------');
                fprintf(1,'\n\n WIE+VIE system solved at freq %.3f MHz, overall time %.2f [sec]\n', ff/1e6, toc(tini));
                fprintf(1,'\n ---------------------------------------------------------------------');
                fprintf(1,'\n ---------------------------------------------------------------------\n\n');

            end

        % -----------------------------------------------------------------
        case 'S' % surface coil model
            
            Nvars = max(COIL.index);
            Nports = length(COIL.port);
            
            % allocate space
            if flags(1)
                ZP = zeros(Nports,Nports,Nfreqs);
            end
            if flags(2)
                Jc = zeros(Nvars,Nports,Nfreqs);
            end
            if flags(3)
                Jb = zeros(L,M,N,3,Nports,Nfreqs);
            end
            if flags(4)
                Sb = zeros(L,M,N,Nports,Nfreqs);
            end
            if flags(5)
                Eb = zeros(L,M,N,3,Nports,Nfreqs);
            end
            if flags(6)
                Bb = zeros(L,M,N,3,Nports,Nfreqs);
            end
            if flags(7)
                Gsar = zeros(Nports,Nfreqs);
            end
            if flags(8)
                Pabs = zeros(Nports,Nfreqs);
            end
            
            for ii = 1:length(freq) % loop on the frequencies
                
                ff = freq(ii);
                
                tini = tic;
                fprintf(1,'\n\n ---------------------------------------------------------------------');
                fprintf(1,'\n ---------------------------------------------------------------------');
                fprintf(1,'\n Solving the SIE+VIE system at freq %.3f MHz\n ', ff/1e6);

                coup = 5; gpu_flag = 1;
                factorized = 1; tolfactor = 1e-5; blockSize = 50;
                % aca_order = 1; aca_tol = 1e-5;
                [Zparam,Jcoil,Jbody] = VSIE_Solver(RHBM,COIL,ff,tol,coup,gpu_flag,factorized,tolfactor,blockSize);

                % store solution
                if flags(1)
                    ZP(:,:,ii) = Zparam;
                end
                if flags(2)
                    Jc(:,:,ii) = Jcoil;
                end
                if sum(flags(3:8))
                    
                    
                    if (flags(5) && flags(9))% need electric field in air
                        % -------------------------------------------------
                        %   Generate circulants if needed
                        if (isempty(RHBM.fN))
                            % we need to compute a new circulant
                            [RHBM.fN] = getOPERATORS(RHBM.r,ff,'N');
                            RHBM.freqfN = ff;
                        else
                            if (RHBM.freqfN ~= ff)
                                [RHBM.fN] = getOPERATORS(RHBM.r,ff,'N');
                                RHBM.freqfN = ff;
                            end
                        end
                    end
                    
                    if flags(6)
                        % -------------------------------------------------
                        %                  Generate circulants if needed
                        if (isempty(RHBM.fK))
                            % we need to compute a new circulant
                            [RHBM.fK] = getOPERATORS(RHBM.r,ff,'K');
                            RHBM.freqfK = ff;
                        else
                            if (RHBM.freqfK ~= ff)
                                [RHBM.fK] = getOPERATORS(RHBM.r,ff,'K');
                                RHBM.freqfK = ff;
                            end
                        end
                    end
                                                                

                    for jj = 1:Nports
                        
                        fprintf(1,'\n\n ----------------');
                        fprintf(1,'\n Port # %d', jj);
                        
                        % compute fields for each channel                    
                        [J,E,B,S,G,P] = VSIE_Fields(RHBM,COIL,Jcoil(:,jj),Jbody(:,jj),ff,flags);

                        if flags(3)
                            Jb(:,:,:,:,jj,ii) = J;
                        end
                        if flags(4)
                            Sb(:,:,:,jj,ii) = S;
                        end
                        if flags(5)
                            Eb(:,:,:,:,jj,ii) = E;
                        end
                        if flags(6)
                            Bb(:,:,:,:,jj,ii) = B;
                        end
                        if flags(7)
                            Gsar(jj,ii) = G;
                        end
                        if flags(8)
                            Pabs(jj,ii) = P;
                        end

                    end

                end
                
                fprintf(1,'\n\n ---------------------------------------------------------------------');
                fprintf(1,'\n\n SIE+VIE system solved at freq %.3f MHz, overall time %.2f [sec]\n', ff/1e6, toc(tini));
                fprintf(1,'\n ---------------------------------------------------------------------');
                fprintf(1,'\n ---------------------------------------------------------------------\n\n');
                
            end
            
        % -----------------------------------------------------------------
        otherwise % wrong coil model or no coil model
            return;
    end

    
    
    
% -------------------------------------------------------------------------
else % there is no human body

    switch COIL.type
        
        case 'W' % wire coil model
        
            Nports = length(COIL.port); % get number of ports
            Nvars = size(COIL.Pcoil,1); % number of segments
            
            % allocate space
            if flags(1)
                ZP = zeros(Nports,Nports,Nfreqs);
            end
            if flags(2)
                Jc = zeros(Nvars,Nports,Nfreqs);
            end
            if flags(5)
                Eb = zeros(L,M,N,3,Nports,Nfreqs);
            end
            if flags(6)
                Bb = zeros(L,M,N,3,Nports,Nfreqs);
            end
            
            for ii = 1:length(freq) % loop on the frequencies
                
                ff = freq(ii);
                
                tini = tic;
                fprintf(1,'\n\n ---------------------------------------------------------------------');
                fprintf(1,'\n ---------------------------------------------------------------------');
                fprintf(1,'\n Solving the WIE system at freq %.3f MHZ\n ', ff/1e6);
                
                [Zparam,Jcoil] = WIE_Solver(COIL,ff);

                if flags(1)
                    ZP(:,:,ii) = Zparam;
                end
                if flags(2)
                    Jc(:,:,ii) = Jcoil;
                end
                if flags(5) || flags(6) % we need to compute the E or B fields
                    for jj = 1:Nports
                        [E,H] = WIE_Radiate(COIL,Jcoil(:,jj),ff,RHBM.r);
                        if flags(5)
                            Eb(:,:,:,:,jj,ii) = E;
                        end
                        if flags(6)
                            Bb(:,:,:,:,jj,ii) = 4*pi*1e-7*H;
                        end
                    end
                end
                
                fprintf(1,'\n\n ---------------------------------------------------------------------');
                fprintf(1,'\n\n WIE system solved at freq %.3f MHz, overall time %.2f [sec]\n', ff/1e6, toc(tini));
                fprintf(1,'\n ---------------------------------------------------------------------');
                fprintf(1,'\n ---------------------------------------------------------------------\n\n');
                
            end

        
        case 'S' % surface coil model
            
            Nvars = max(COIL.index);
            Nports = length(COIL.port);
            
            % allocate space
            if flags(1)
                ZP = zeros(Nports,Nports,Nfreqs);
            end
            if flags(2)
                Jc = zeros(Nvars,Nports,Nfreqs);
            end
            if flags(5)
                Eb = zeros(L,M,N,3,Nports,Nfreqs);
            end
            if flags(6)
                Bb = zeros(L,M,N,3,Nports,Nfreqs);
            end
            
            for ii = 1:length(freq) % loop on the frequencies
                
                ff = freq(ii);
                
                tini = tic;
                fprintf(1,'\n\n ---------------------------------------------------------------------');
                fprintf(1,'\n ---------------------------------------------------------------------');
                fprintf(1,'\n Solving the SIE system at freq %.3f MHz\n ', ff/1e6);
                
                
                [Zparam,Jcoil] = SIE_Solver(COIL,ff);
                if flags(1)
                    ZP(:,:,ii) = Zparam;
                end
                if flags(2)
                    Jc(:,:,ii) = Jcoil;
                end
                if flags(5) || flags(6) % we need to compute the E or B fields
                    for jj = 1:Nports
                        [E,H] = SIE_Radiate(COIL,Jcoil(:,jj),ff,RHBM.r);
                        if flags(5)
                            Eb(:,:,:,:,jj,ii) = E;
                        end
                        if flags(6)
                            Bb(:,:,:,:,jj,ii) = 4*pi*1e-7*H;
                        end
                    end
                end
                
                fprintf(1,'\n\n ---------------------------------------------------------------------');
                fprintf(1,'\n\n SIE system solved at freq %.3f MHz, overall time %.2f [sec]\n', ff/1e6, toc(tini));
                fprintf(1,'\n ---------------------------------------------------------------------');
                fprintf(1,'\n ---------------------------------------------------------------------\n\n');
                
            end
            
            
        otherwise % wrong coil model or no coil model
            return;            
    end

end

