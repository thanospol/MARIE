function [epsilon_r,sigma_e,rho,mu_r,sigma_m] = RHBM_Load(epsilonfile,sigmaefile,densityfile,mufile,sigmamfile,nX,nY,nZ)
%%    Function to load material properties
% _________________________________________________________________________
%
%       Loads data from voxelized files
%
% _________________________________________________________________________
%
%% INPUT
%           files where the density, epsilon and sigma are stored
%           nX,nY,nZ    number of elements in each direction
%           f           frequency
%
%% OUTPUT:
%           properties
%
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
%                              Loda file data
% -------------------------------------------------------------------------


if ~isempty(epsilonfile)
    file=fopen(epsilonfile,'rb');
    if (file ~= -1)
        epsilon_r = reshape((fread(file,'float32')),nX,nY,nZ);
        fclose(file);
    else
        epsilon_r = ones(nX,nY,nZ);
    end
else
    epsilon_r = ones(nX,nY,nZ);
end

if ~isempty(sigmaefile)
    file=fopen(sigmaefile,'rb');
    if (file ~= -1)
        sigma_e = reshape((fread(file,'float32')),nX,nY,nZ);
        fclose(file);
    else
        sigma_e = zeros(nX,nY,nZ);
    end
else
    sigma_e = zeros(nX,nY,nZ);
end

if ~isempty(densityfile)
    file=fopen(densityfile,'rb');
    if (file ~= -1)
        rho = reshape((fread(file,'float32')),nX,nY,nZ);
        fclose(file);
    else
        rho = [];
    end
else
    rho = [];
end

if ~isempty(mufile)
    file=fopen(mufile,'rb');
    if (file ~= -1)
        mu_r = reshape((fread(file,'float32')),nX,nY,nZ);
        fclose(file);
    else
        mu_r = [];
    end
else
    mu_r = [];
end


if ~isempty(sigmamfile)
    file=fopen(sigmamfile,'rb');
    if (file ~= -1)
        sigma_m = reshape((fread(file,'float32')),nX,nY,nZ);
        fclose(file);
    else
        sigma_m = [];
    end
else
    sigma_m = [];
end

clear file

% % % -------------------------------------------------------------------------
% % %                             Remove shitty stuff
% % % -------------------------------------------------------------------------
% % 
% % for mx = 1:nX
% %     for my = 1:nY
% %         for mz = 1:nZ
% %             if sigma_e(mx,my,mz) > 140
% %                 sigma_e(mx,my,mz) = 0.0;
% %                 epsilon_r(mx,my,mz) = 1.0;
% %                 rho(mx,my,mz) = 0.0;
% %             end
% %         end
% %     end
% % end




