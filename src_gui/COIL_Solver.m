function [ZP,Jc] = COIL_Solver(COIL,freq)
%%    Simple driver to solve a frequency sweep for a coil model
% _________________________________________________________________________
%
%       Solves the frequency domain analysis of a coil model
%
% _________________________________________________________________________
%
%% INPUT
%   COIL       Coil model structure
%   freq       Frequency vector (in Hz)
%
%
%% OUTPUT
%   ZP          Z parameter matrix (NpxNpxNfreqs)
%   Jc          Coil basis (current) coefficients (NcxNpxNfreqs)
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

ZP = [];
Jc = [];

if(nargin < 2 )
   fprintf(1, '\n ERROR: not enough arguments\n');
   return
end

% -------------------------------------------------------------------------
%   check the existence of different cases and solve
% -------------------------------------------------------------------------

freq = squeeze(freq);
Nfreqs = length(freq); % number of points


switch COIL.type
    
    case 'W' % wire coil model
        
        Nports = length(COIL.port); % get number of ports
        Nvars = size(COIL.Pcoil,1); % number of segments
        
        tini = tic;
        fprintf(1,'\n\n ---------------------------------------------------------------------');
        fprintf(1,'\n ---------------------------------------------------------------------');
        fprintf(1,'\n WIE system frequency sweep\n ');
        
        % allocate space
        ZP = zeros(Nports,Nports,Nfreqs);
        Jc = zeros(Nvars,Nports,Nfreqs);
        
        for ii = 1:length(freq) % loop on the frequencies
            
            ff = freq(ii);
            [Zparam,Jcoil] = WIE_Solver(COIL,ff);
            
            ZP(:,:,ii) = Zparam;
            Jc(:,:,ii) = Jcoil;
            
        end
        
        fprintf(1,'\n\n ---------------------------------------------------------------------');
        fprintf(1,'\n\n WIE frequency sweep for %d freqs. done, overall time %.2f [sec]\n', ii, toc(tini));
        fprintf(1,'\n ---------------------------------------------------------------------');
        fprintf(1,'\n ---------------------------------------------------------------------\n\n');
        
        
    case 'S' % surface coil model
        
        Nvars = max(COIL.index);
        Nports = length(COIL.port);
        
        tini = tic;
        fprintf(1,'\n\n ---------------------------------------------------------------------');
        fprintf(1,'\n ---------------------------------------------------------------------');
        fprintf(1,'\n SIE system frequency sweep\n ');
        
        % allocate space
        ZP = zeros(Nports,Nports,Nfreqs);
        Jc = zeros(Nvars,Nports,Nfreqs);
        
        for ii = 1:length(freq) % loop on the frequencies
            
            ff = freq(ii);
            [Zparam,Jcoil] = SIE_Solver(COIL,ff);
            
            ZP(:,:,ii) = Zparam;
            Jc(:,:,ii) = Jcoil;
            
        end
        
        fprintf(1,'\n\n ---------------------------------------------------------------------');
        fprintf(1,'\n\n SIE frequency sweep for %d freqs. done, overall time %.2f [sec]\n', ii, toc(tini));
        fprintf(1,'\n ---------------------------------------------------------------------');
        fprintf(1,'\n ---------------------------------------------------------------------\n\n');
        
        
    otherwise % wrong coil model or no coil model
        return;
end



