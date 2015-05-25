function [x,y,z,epsilon_r,sigma_e,rho,mu_r,sigma_m,material] = RHBMDat_Parse(voxelfile,xfile,yfile,zfile)
%%    Function to load material properties
% _________________________________________________________________________
%
%       Reads file with voxel properties and coordinates
%       Returns the propertie
%
% _________________________________________________________________________
%
%% INPUT
%   voxefile    File with LxMxN entries with the properties
%               Format of the file is
%               xindex yindex zindex material_number epsilon_r sigma_e mu_r sigma_h rho
%   xfile       File with L values of the x coordinates for the 3D grid
%   yfile       File with M values of the y coordinates for the 3D grid
%   zfile       File with N values of the z coordinates for the 3D grid
%
%
%% OUTPUT
%   x           x coordinates for the 3D grid (Lx1)
%   y           y coordinates for the 3D grid (Mx1)
%   z           z coordinates for the 3D grid (Nx1)
%   epsilon_r   relative epsilon (LxMxN)
%   sigma_e     electric  (LxMxN)
%   rho         density (LxMxN)
%   mu_r        relative mu (LxMxN)
%   sigma_m     magnetic (LxMxN)
%   material    number of the corresponding material (LxMxN)
%
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
% Read coordinate files
% -------------------------------------------------------------------------

tini = tic;
fprintf(1, '\n Loading "%s"  with\n   %s\n   %s\n   %s\n', voxelfile,xfile,yfile,zfile);
fprintf(1, ' as coordinate files... this operation may take some minutes\n');

% allocate space for coordinates
x = zeros(1000,1);
y = zeros(1000,1);
z = zeros(1000,1);

% initialize counters
xcount = 0;
ycount = 0;
zcount = 0;


% open file
fx = fopen(xfile, 'r');
% loop to read files
while ~(feof(fx)) % not finished with file fx
       
    xline = fgetl(fx); % read line
    
    if ~(xline(1) == '%') % not a comment
        xcount = xcount + 1;
        x(xcount) = sscanf(xline,'%f'); % read coordinate
    end
    
end
fclose(fx);    
fprintf(1, '.');

% open file
fy = fopen(yfile, 'r');
% loop to read files
while ~(feof(fy)) % not finished with fy
    
    yline = fgetl(fy); % read line
    
    if (yline(1) ~= '%') % not a comment
        ycount = ycount + 1;
        y(ycount) = sscanf(yline,'%f'); % read coordinate
    end
    
end
fclose(fy);
fprintf(1, '.');

% open file
fz = fopen(zfile, 'r');
% loop to read files
while ~(feof(fz)) % not finished with
    
    zline = fgetl(fz); % read line
    
    if ~(zline(1) == '%') % not a comment
        zcount = zcount + 1;
        z(zcount) = sscanf(zline,'%f'); % read coordinate
    end
    
end
fclose(fz);
fprintf(1, '.');

% cut to final size
x = x(1:xcount);
y = y(1:ycount);
z = z(1:zcount);

% number of voxels
L = length(x);
M = length(y);
N = length(z);

% voxelized
fvox = fopen(voxelfile, 'r');

% allocate space
material = zeros(L,M,N);
epsilon_r = zeros(L,M,N);
sigma_e = zeros(L,M,N);
mu_r = zeros(L,M,N);
sigma_m	= zeros(L,M,N);
rho = zeros(L,M,N);

count = 0;
while ~(feof(fvox)) % not finished with
        
        line = fgetl(fvox); % read line
        
        if (line(1) ~= '%') % not a comment
            count = count + 1;

            vec = sscanf(line,'%f'); % read coordinate

            xidx = vec(1)+1; % x index
            yidx = vec(2)+1; % y index
            zidx = vec(3)+1; % z index
            
            material(xidx,yidx,zidx) = vec(4); % material number
            
            % properties
            epsilon_r(xidx,yidx,zidx) = vec(5);
            sigma_e(xidx,yidx,zidx) = vec(6);
            mu_r(xidx,yidx,zidx) = vec(7);
            sigma_m(xidx,yidx,zidx) = vec(8);
            rho(xidx,yidx,zidx) = vec(9);            
            
        end
        
        if (rem(count,500000) == 0)
            fprintf(1, '.');
        end

end

fprintf(1, '\n');
fclose(fvox);

if (count ~= L*M*N)
    fprintf(1, '\n\n Warning: number of voxels loaded and number of coordinates do not match\n\n');
end

fprintf(1, '\n Load of %s  done\n %d lines read, %g seconds\n', voxelfile,count, toc(tini));

