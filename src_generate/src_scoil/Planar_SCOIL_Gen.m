function [ok] = Planar_SCOIL_Gen(filename,Ncoils,Tports,Rports,Bports,Lports,Rad,Xpos,Ypos,Zpos,Len,Ssp,Wid,Gap,IniAngle,ShieldRad,ShieldLen)
%%    Surface Planar Coil Array Generator
% _________________________________________________________________________
%
%
%   Def.: Fucntion 
%
% _________________________________________________________________________
%
%
%% INPUT
%       filename - name of the geometry file to load data
%       Ncoils - number of coils in array
%       Tports - array with position of top ports
%       Rports - array with position of right side ports
%       Bports - array with position of bottom ports
%       Lports - array with position of left side ports
%       Rad - radius: the distance of the coil array to the center
%       Xpos,Ypos,Zpos - position of the center of coil array
%       Ssp - external span of each coil
%       Len - external length in z direction of each coil
%       Wid - width of the trace of the coil
%       Gap - gap of the ports
%       IniAngle - initial angle for the center of the first coil
%
%
%% OUTPUT
%       ok
%
%
% -------------------------------------------------------------------------
%
%   J. Fernandez Villena -- jvillena@mit.edu
%   A.G. Polimeridis -- thanos_p@mit.edu
%   Computational Prototyping Group, RLE at MIT
%
% _________________________________________________________________________
 

% check values
if (Gap <= 0)
    Gap = 1e-6;
end

if (Wid <= 0)
    Wid = 5e-3;
end

if (Len <= 0) 
    Len = 0.12;
end
  
if (Ssp <= 0)
    Ssp = 0.08;
end


% generate the name for the geometry file
geofilename = sprintf('%s.geo', filename);
fid = fopen(geofilename, 'w');
fprintf(fid, '\n\n\n // --------------------------------------------------------------- '); 
fprintf(fid, '\n // --------------------------------------------------------------- ');
fprintf(fid, '\n // Planar COIL ARRAY with %d COILS', Ncoils); 
fprintf(fid, '\n // geometry created using MARIE coil generator');  
fprintf(fid, '\n // --------------------------------------------------------------- ');
fprintf(fid, '\n // --------------------------------------------------------------- '); 
fclose(fid);


% loop on coils and call the coil generator
DeltaAngle = 360/Ncoils;
for ii = 1:Ncoils
    
    Arot = ((ii-1)*DeltaAngle+IniAngle)*pi/180; % rotation angle
    [ok] = gen_coil_planarpatches(geofilename,ii,Tports,Rports,Bports,Lports,Rad,Xpos,Ypos,Zpos,Arot,Len,Ssp,Wid,Gap);
    
end
    

% add Shield if exist
if (~isempty(ShieldRad)) && (~isempty(ShieldLen))
    [ok] = Shield_SCOIL_Gen(geofilename,ShieldRad,Xpos,Ypos,Zpos,ShieldLen);
end


% call the gmsh
if ispc
    command = sprintf('.\\src_generate\\src_scoil\\gmsh %s', geofilename);
else
    if ismac
        command = sprintf('./src_generate/src_scoil/gmsh.app/Contents/MacOS/gmsh %s', geofilename);
    end
end
[ok, result] = system(command);

