function [PULSE] = pulse_apply(SOL,COIL,RHBM,pulsefile)
%%    Simple driver to apply pulse to a solution
% _________________________________________________________________________
%
%       Applies the pulse from a .mps file (marie pulse sequence) to a SOL
%       returns the structure with the pulse solution
%
% _________________________________________________________________________
%
%% INPUT
%   SOL        Solution structure
%   COIL       Coil model structure
%   RHBM       Body model structure
%   pulsefile  name of the pulse file
%
%
%% OUTPUT
%   PULSE      structure with the pulse solution
%
% -------------------------------------------------------------------------
%
%%   This function is part of MARIE
%   MARIE - Magnetic Resonance Integral Equation suite
%           Jorge Fernandez Villena   -- jvillena@mit.edu
%           Athanasios G. Polimeridis -- thanos_p@mit.edu
%           Copyright © 2014
%           RLE Computational Prototyping Group, MIT
% 
%           This software is free and open source
%           Distributed under the GNU-GPLv3 terms
%           For details see MARIE_license.txt
%
% _________________________________________________________________________



% -------------------------------------------------------------------------
% Initialization of variables
% -------------------------------------------------------------------------

if(nargin < 3 )
   fprintf(1, '\n ERROR: not enough arguments\n');
   return
end


tini = tic;
fprintf(1,'\n\n ---------------------------------------------------------------------');
fprintf(1,'\n ---------------------------------------------------------------------');
fprintf(1,'\n Applying pulse sequence from:\n  %s\n ', pulsefile);

% -------------------------------------------------------------------------
%   parse the mps file and get excitation data
% -------------------------------------------------------------------------

[Npulses,Mports,Spulse,Vpulse,Ipulse] = import_mps(pulsefile);

% -------------------------------------------------------------------------
%   preallocate the variables of PULSE
% -------------------------------------------------------------------------
% 

[Nvars,Nports,Nfreqs] = size(SOL.Jcoil); % number of coil variables
[L,M,N,~,~,~] = size(SOL.Jsol); % size of the domain

% check consistency
if (Nports ~= Mports)
    PULSE = [];
    fprintf(1, '\n ERROR: number of ports does not match\n');
    return;
end

PULSE.Jcoil = zeros(Nvars,Npulses,Nfreqs);
PULSE.Jsol = zeros(L,M,N,3,Npulses,Nfreqs);
PULSE.Ssol = zeros(L,M,N,Npulses,Nfreqs);
PULSE.Esol = zeros(L,M,N,3,Npulses,Nfreqs);
PULSE.Bsol = zeros(L,M,N,3,Npulses,Nfreqs);
PULSE.Pabs = zeros(Npulses,Nfreqs);
PULSE.Gsar = zeros(Npulses,Nfreqs);
PULSE.step = Spulse;
PULSE.freq = SOL.freq;
PULSE.Ipulse = Ipulse;
PULSE.Vpulse = Vpulse;

% -------------------------------------------------------------------------
%   define the kind of excitation depending on the coil model
% -------------------------------------------------------------------------

switch (COIL.type)
    case 'S'
        % for surface coils, the solution is computed for unitary currents
        % so we need to scale by the currents at each port
        WP = Ipulse;
    case 'W'
        % for wire coils, the solution is computed for unitary voltages
        % so we need to scale by the voltages at each port
        WP = Vpulse;
    otherwise
        % no coil, error
        PULSE = [];
        fprintf(1, '\n ERROR: there is no available coil model\n');
        return;
end

% -------------------------------------------------------------------------
%   Apply the corresponding pulse weighting
% -------------------------------------------------------------------------

% dimensions
dx = RHBM.r(2,1,1,1) - RHBM.r(1,1,1,1);
Gram = dx^3;

% loop on the frequencies for weighting the data
for jj = 1:Nfreqs
    
    % get the solution to each pulse for each frequency
    for ii = 1:Npulses
        
        for kk = 1:Nports
            PULSE.Jcoil(:,ii,jj) = PULSE.Jcoil(:,ii,jj)+SOL.Jcoil(:,kk,jj)*WP(kk,ii);
            PULSE.Jsol(:,:,:,:,ii,jj) = PULSE.Jsol(:,:,:,:,ii,jj)+SOL.Jsol(:,:,:,:,kk,jj)*WP(kk,ii);
            PULSE.Esol(:,:,:,:,ii,jj) = PULSE.Esol(:,:,:,:,ii,jj)+SOL.Esol(:,:,:,:,kk,jj)*WP(kk,ii);
            PULSE.Bsol(:,:,:,:,ii,jj) = PULSE.Bsol(:,:,:,:,ii,jj)+SOL.Bsol(:,:,:,:,kk,jj)*WP(kk,ii);
            PULSE.Ssol(:,:,:,ii,jj) = PULSE.Ssol(:,:,:,ii,jj)+SOL.Ssol(:,:,:,kk,jj)*WP(kk,ii);        
        end
        
        % Pabs and GSAR
        Pabs = PULSE.Esol(:,:,:,:,ii,jj).*conj(PULSE.Jsol(:,:,:,:,ii,jj));
        Pabs = sum(Pabs,4);
        Pabs = 0.5*real(sum(Pabs(:)))*Gram;
        PULSE.Pabs(ii,jj) = Pabs;
        PULSE.Gsar(ii,jj) = sum(PULSE.Ssol(:))*Gram;          
    end  
end

% -------------------------------------------------------------------------
%   and done
% -------------------------------------------------------------------------


fprintf(1,'\n\n Pulse sequence applied, overall time %.2f [sec]', toc(tini));
fprintf(1,'\n ---------------------------------------------------------------------\n\n');


