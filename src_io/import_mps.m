function [Npulses,Mports,Spulse,Vpulse,Ipulse] = import_mps(filename)
%%    Import voltage and current from a marie pulse sequence file
% _________________________________________________________________________
%
%
%   Def.: Loads the voltages and currents at ports due to sequence
%
% _________________________________________________________________________
%
%
%% INPUT
%       filename - name of the mps file with the data
%
%
%% OUTPUT
%       Npulses - number of pulses applied
%       Mports - number of ports
%       Spulse - Npulses array with the pulse step or time
%       Vpulse - MportsxNpulses with voltage at all ports for each pulse
%       Ipulse - MportsxNpulses with currents at all ports for each pulse
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

if(nargin < 1 )
   fprintf(1, '\n ERROR: not enough arguments\n');
   return
end


% -------------------------------------------------------------------------
% Open file, read data and store parameters
% -------------------------------------------------------------------------

fid = fopen(filename, 'r'); % check of mps file is done outside


% parse the comments (start with !, until find the data, starting with #)
while ~(feof(fid)) % not finished with
        
        line = fgetl(fid); % read line
        
        if (line(1) == '#') % mark to start of data
           
            vec = sscanf(line(2:end),'%f'); % read number of pulses and of ports

            Npulses = vec(1); % number of pulses to apply
            Mports = vec(2); % number of ports of the data
            
            break;
            
        end
        
end

% allocate
Spulse = zeros(Npulses,1);
Vpulse = zeros(Mports,Npulses);
Ipulse = zeros(Mports,Npulses);

% now read the data
for ii = 1:Npulses
    
    %read first line that has pulse number or time step
    line = fgetl(fid); % read line
    vec = sscanf(line,'%f'); % read pulse step and port excitation
                
    Spulse(ii) =  vec(1);
    Vpulse(1,ii) = vec(2) + 1j*vec(3);
    Ipulse(1,ii) = vec(4) + 1j*vec(5);
    
    for jj = 2:Mports % remaining ports
        
        line = fgetl(fid); % read line
        vec = sscanf(line,'%f'); % read port excitation
        
        Vpulse(jj,ii) = vec(1) + 1j*vec(2);
        Ipulse(jj,ii) = vec(3) + 1j*vec(4);
        
    end
    
end
        
      
% -------------------------------------------------------------------------
% and it is done
% ------------------------------------------------------------------------- 
        
fclose(fid);



