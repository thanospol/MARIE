function [ok] = Birdcage_SCOIL_Gen(filename,Nleg,Rad,Xpos,Ypos,Zpos,Len,Wid,Gap,IniAngle,ShieldRad,ShieldLen)
%%    Surface Conformal Birdcage Generator
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
%       Nleg - number of legs
%       Rad - radius: the distance of the coil array to the center
%       Xpos,Ypos,Zpos - position of the center of coil array
%       Len - external length in z direction of each coil
%       Wid - width of the trace of the coil
%       Gap - gap of the ports
%       IniAngle - initial angle for the center of the first port
%
%
%% OUTPUT
%       ok
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
 

% check values
if (Gap <= 0)
    Gap = 1e-6;
end

if (Wid <= 0)
    Wid = 5e-3;
end

if (Len <= 0) 
    Len = 0.18;
end
  

% generate the name for the geometry file
geofilename = sprintf('%s.geo', filename);
fid = fopen(geofilename, 'w');
fprintf(fid, '\n\n\n // --------------------------------------------------------------- '); 
fprintf(fid, '\n // --------------------------------------------------------------- ');
fprintf(fid, '\n // BIRDCAGE with %d legs', Nleg); 
fprintf(fid, '\n // geometry created using MARIE coil generator');  
fprintf(fid, '\n // --------------------------------------------------------------- ');
fprintf(fid, '\n // --------------------------------------------------------------- '); 
fclose(fid);

Res = 3*Wid;
Arot = IniAngle*pi/180; % rotation of the coil

% -------------------------------------------------------------------------
% generate the data
% -------------------------------------------------------------------------

Asp = 2*pi/Nleg; % external angle span in radians
len = Len - 2*Wid; % internal length
wsp = asin(Wid/Rad); % width span in radians
asp = Asp - 2*wsp; % internal angle span in rads
gsp = asin(Gap/Rad);% gap span in rads


% local resolution
Rescorner = Res;
Respatch = Res;
Resport = Res;

% we will have N structures:
% point array
% line structure: one for each line or circle, 

pcount = 0;
lcount = 0;

localexternal = [ 2; 4; 6; 8; 10; 12; 14; 16; 17; 20; 21; 22];
external = [];

Points = [];
Lines = [];
Patch = [];

% -------------------------------------------------------------------------
% Center points for curves

x = 0; y = 0; z = +Len/2; % center for top curve
Points = [Points; x, y, z, Res];
pcount = pcount + 1;
x = 0; y = 0; z = +len/2; % center for top internal top curve
Points = [Points; x, y, z, Res];
pcount = pcount + 1;
x = 0; y = 0; z = -len/2; % center for internal bottom curve
Points = [Points; x, y, z, Res];
pcount = pcount + 1;
x = 0; y = 0; z = -Len/2; % center for top external bottom curve
Points = [Points; x, y, z, Res];
pcount = pcount + 1;




% loop to generate the different elements of the birdcage coil

for legcount = 1:Nleg
    
    
    inialpha = Asp*(legcount-1);
    
    % point distribution in birdcage element is going to be
    %
    %           2       3   6      7
    %   p1+                             p2-
    %           1       4   5      8
    %
    %
    %
    %           14      15  10      11
    %   pN+1+                           pN+2-
    %           13       16  9      12
    %
    
    % -------------------------------------------------------------------------
    % Generate the points
    
    pshift = pcount;
    lshift = lcount;
    
    alpha = inialpha+gsp/2;
    x = Rad*sin(alpha); y = Rad*cos(alpha); z = +len/2; % 1
    Points = [Points; x, y, z, Res];
    pcount = pcount + 1;
    x = Rad*sin(alpha); y = Rad*cos(alpha); z = +Len/2; % 2
    Points = [Points; x, y, z, Res];
    pcount = pcount + 1;
    
    alpha = inialpha+(Asp-wsp)/2;
    x = Rad*sin(alpha); y = Rad*cos(alpha); z = +Len/2; % 3
    Points = [Points; x, y, z, Res];
    pcount = pcount + 1;
    x = Rad*sin(alpha); y = Rad*cos(alpha); z = +len/2; % 4
    Points = [Points; x, y, z, Res];
    pcount = pcount + 1;
    
    alpha = inialpha+(Asp+wsp)/2;
    x = Rad*sin(alpha); y = Rad*cos(alpha); z = +len/2; % 5
    Points = [Points; x, y, z, Res];
    pcount = pcount + 1;
    x = Rad*sin(alpha); y = Rad*cos(alpha); z = +Len/2; % 6
    Points = [Points; x, y, z, Res];
    pcount = pcount + 1;
    
    alpha = inialpha+Asp-gsp/2;
    x = Rad*sin(alpha); y = Rad*cos(alpha); z = +Len/2; % 7
    Points = [Points; x, y, z, Res];
    pcount = pcount + 1;
    x = Rad*sin(alpha); y = Rad*cos(alpha); z = +len/2; % 8
    Points = [Points; x, y, z, Res];
    pcount = pcount + 1;
    
    alpha = inialpha+(Asp+wsp)/2;
    x = Rad*sin(alpha); y = Rad*cos(alpha); z = -Len/2; % 9
    Points = [Points; x, y, z, Res];
    pcount = pcount + 1;
    x = Rad*sin(alpha); y = Rad*cos(alpha); z = -len/2; % 10
    Points = [Points; x, y, z, Res];
    pcount = pcount + 1;
    
    alpha = inialpha+Asp-gsp/2;
    x = Rad*sin(alpha); y = Rad*cos(alpha); z = -len/2; % 11
    Points = [Points; x, y, z, Res];
    pcount = pcount + 1;
    x = Rad*sin(alpha); y = Rad*cos(alpha); z = -Len/2; % 12
    Points = [Points; x, y, z, Res];
    pcount = pcount + 1;
    
    alpha = inialpha+gsp/2;
    x = Rad*sin(alpha); y = Rad*cos(alpha); z = -Len/2; % 13
    Points = [Points; x, y, z, Res];
    pcount = pcount + 1;
    x = Rad*sin(alpha); y = Rad*cos(alpha); z = -len/2; % 14
    Points = [Points; x, y, z, Res];
    pcount = pcount + 1;
    
    alpha = inialpha+(Asp-wsp)/2;
    x = Rad*sin(alpha); y = Rad*cos(alpha); z = -len/2; % 15
    Points = [Points; x, y, z, Res];
    pcount = pcount + 1;
    x = Rad*sin(alpha); y = Rad*cos(alpha); z = -Len/2; % 16
    Points = [Points; x, y, z, Res];
    pcount = pcount + 1;
    
    
    % -------------------------------------------------------------------------
    % Top Left patch
    
    idxp1 = pshift+1;
    idxp2 = pshift+2;
    idxp3 = pshift+3;
    idxp4 = pshift+4;
    
    % lines --------------------------------------------------
    Lines = [Lines; idxp1, idxp2, -100];
    lcount = lcount + 1;
    idxl1 = lcount;
    Lines = [Lines; idxp2, idxp3, 1];
    lcount = lcount + 1;
    idxl2 = lcount;
    Lines = [Lines; idxp3, idxp4, -100];
    lcount = lcount + 1;
    idxl3 = lcount;
    Lines = [Lines; idxp4, idxp1, 2];
    lcount = lcount + 1;
    idxl4 = lcount;
    
    % surface --------------------------------------------------
    Patch = [Patch; idxl1 idxl2 idxl3 idxl4];
    
    
    % -------------------------------------------------------------------------
    % Top Right patch
    
    % points --------------------------------------------------
    idxp1 = pshift+5;
    idxp2 = pshift+6;
    idxp3 = pshift+7;
    idxp4 = pshift+8;
    
    % lines --------------------------------------------------
    Lines = [Lines; idxp1, idxp2, -100];
    lcount = lcount + 1;
    idxl1 = lcount;
    Lines = [Lines; idxp2, idxp3, 1];
    lcount = lcount + 1;
    idxl2 = lcount;
    Lines = [Lines; idxp3, idxp4, -100];
    lcount = lcount + 1;
    idxl3 = lcount;
    Lines = [Lines; idxp4, idxp1, 2];
    lcount = lcount + 1;
    idxl4 = lcount;
    
    % surface --------------------------------------------------
    Patch = [Patch; idxl1 idxl2 idxl3 idxl4];
    
    
    % -------------------------------------------------------------------------
    % Bottom Right  patch
    
    % points --------------------------------------------------
    idxp1 = pshift+9;
    idxp2 = pshift+10;
    idxp3 = pshift+11;
    idxp4 = pshift+12;
    
    % lines --------------------------------------------------
    Lines = [Lines; idxp1, idxp2, -100];
    lcount = lcount + 1;
    idxl1 = lcount;
    Lines = [Lines; idxp2, idxp3, 3];
    lcount = lcount + 1;
    idxl2 = lcount;
    Lines = [Lines; idxp3, idxp4, -100];
    lcount = lcount + 1;
    idxl3 = lcount;
    Lines = [Lines; idxp4, idxp1, 4];
    lcount = lcount + 1;
    idxl4 = lcount;
    
    % surface --------------------------------------------------
    Patch = [Patch; idxl1 idxl2 idxl3 idxl4];
    
    
    % -------------------------------------------------------------------------
    % Bottom Left patch
    
    % points --------------------------------------------------
    idxp1 = pshift+13;
    idxp2 = pshift+14;
    idxp3 = pshift+15;
    idxp4 = pshift+16;
    
    % lines --------------------------------------------------
    Lines = [Lines; idxp1, idxp2, -100];
    lcount = lcount + 1;
    idxl1 = lcount;
    Lines = [Lines; idxp2, idxp3, 3];
    lcount = lcount + 1;
    idxl2 = lcount;
    Lines = [Lines; idxp3, idxp4, -100];
    lcount = lcount + 1;
    idxl3 = lcount;
    Lines = [Lines; idxp4, idxp1, 4];
    lcount = lcount + 1;
    idxl4 = lcount;
    
    % surface --------------------------------------------------
    Patch = [Patch; idxl1 idxl2 idxl3 idxl4];
    
    
    % -------------------------------------------------------------------------
    % Top central patch
    
    idxp1 = pshift+4;
    idxp2 = pshift+3;
    idxp3 = pshift+6;
    idxp4 = pshift+5;
    
    % lines --------------------------------------------------
    idxl1 = lshift+3; % line is already there
    Lines = [Lines; idxp2, idxp3, 1];
    lcount = lcount + 1;
    idxl2 = lcount;
    idxl3 = lshift+5; % line is already there
    Lines = [Lines; idxp4, idxp1, 2];
    lcount = lcount + 1;
    idxl4 = lcount;
    
    % surface --------------------------------------------------
    Patch = [Patch; -idxl1 idxl2 -idxl3 idxl4];
    
    
    % -------------------------------------------------------------------------
    % Bottom central patch
    
    idxp1 = pshift+16;
    idxp2 = pshift+15;
    idxp3 = pshift+10;
    idxp4 = pshift+9;
    
    % lines --------------------------------------------------
    idxl1 = lshift+15; % line is already there
    Lines = [Lines; idxp2, idxp3, 3];
    lcount = lcount + 1;
    idxl2 = lcount;
    idxl3 = lshift+9; % line is already there
    Lines = [Lines; idxp4, idxp1, 4];
    lcount = lcount + 1;
    idxl4 = lcount;
    
    % surface --------------------------------------------------
    Patch = [Patch; -idxl1 idxl2 -idxl3 idxl4];
    
    
    % -------------------------------------------------------------------------
    % Leg patch
    
    idxp1 = pshift+15;
    idxp2 = pshift+4;
    idxp3 = pshift+5;
    idxp4 = pshift+10;
    
    % lines --------------------------------------------------
    Lines = [Lines; idxp1, idxp2, -100];
    lcount = lcount + 1;
    idxl1 = lcount;
    idxl2 = lshift+18; % line is already there
    Lines = [Lines; idxp3, idxp4, -100];
    lcount = lcount + 1;
    idxl3 = lcount;
    idxl4 = lshift+19; % line is already there
    
    % surface --------------------------------------------------
    Patch = [Patch; idxl1 -idxl2 idxl3 -idxl4];
    
    
    % define external lines for element
    external = [external; localexternal+lshift];
    
end

% -------------------------------------------------------------------------
% Now define the ports

Ports = zeros(2*Nleg,2);

% Initial leg top port
Ports(1,1) = 1;
Ports(1,2) = 22*(Nleg-1)+7;

% Initial leg bottom port
Ports(Nleg+1,1) = 13;
Ports(Nleg+1,2) = 22*(Nleg-1)+11;


for legcount = 1:Nleg-1
    
    portshift = 22*legcount;
    prevshift = 22*(legcount-1);
    
    % top port
    Ports(legcount+1,1) = portshift+1;
    Ports(legcount+1,2) = prevshift+7;
    
    % bottom port
    Ports(Nleg+legcount+1,1) = portshift+13;
    Ports(Nleg+legcount+1,2) = prevshift+11;
    
end



% Print all the stuff in the file
fid = fopen(geofilename, 'w');

% fprintf(fid, '\n Mesh.Format      = 1; // 1=.msh format by default');
% fprintf(fid, '\n Mesh.Algorithm   = 2; // 2D mesh algorithm (1=MeshAdapt, 5=Delaunay, 6=Frontal) Default value: 1 ');



fprintf(fid, '\n\n\n // --------------------------------------------------------------- '); 
fprintf(fid, '\n // Birdcage with %d Legs', Nleg);  
fprintf(fid, '\n //        Radius       %g m', Rad);
fprintf(fid, '\n //        Rotation     %g rads', Arot);
fprintf(fid, '\n //        Translation Z  %g m', Zpos);
fprintf(fid, '\n //        Translation Y  %g m', Ypos);
fprintf(fid, '\n //        Translation X  %g m', Xpos);

fprintf(fid, '\n\n // Resolution');
fprintf(fid, '\n RES = %f;', Res);

Cnum = 1;
Nstart = Cnum*1000;
Pstart = Nstart + 100;
Lstart = Pstart + 400;
Ostart = Lstart + 200;
Sstart = Ostart + 100;
Tstart = Sstart + 100;

fprintf(fid, '\n\n  ');


% print the points
for ii = 1:pcount

    x = Points(ii,1)*cos(Arot) - Points(ii,2)*sin(Arot) + Xpos;
    y = Points(ii,1)*sin(Arot) + Points(ii,2)*cos(Arot) + Ypos;
    z = Points(ii,3) + Zpos;    
    
    if (abs(x) < eps)
        x = 0;
    end
    if (abs(y) < eps)
        y = 0;
    end
    if (abs(z) < eps)
        z = 0;
    end
    
    fprintf(fid, '\n P%d = %d;',Pstart+ii, Pstart+ii);
    fprintf(fid, '\n Point(P%d) = {%2.16g, %2.16g, %2.16g, %2.16g};',Pstart+ii, x, y, z, Points(ii,4));
    
end


fprintf(fid, '\n\n  ');


% print the lines
for ii = 1:lcount

    if (Lines(ii,3) < 0) % straight line
         
        fprintf(fid, '\n L%d = %d;',Lstart+ii, Lstart+ii);
        fprintf(fid, '\n Line(L%d) = {P%d, P%d};',Lstart+ii, Pstart+Lines(ii,1), Pstart+Lines(ii,2));
    
    else
        
        fprintf(fid, '\n L%d = %d;',Lstart+ii, Lstart+ii);
        fprintf(fid, '\n Circle(L%d) = {P%d, P%d, P%d};',Lstart+ii, Pstart+Lines(ii,1), Pstart+Lines(ii,3), Pstart+Lines(ii,2)); 
        
    end
        
    
end


fprintf(fid, '\n\n  ');
% print the loops 
for ii = 1:size(Patch,1)
    
    fprintf(fid, '\n O%d = %d;', Ostart+ii, Ostart+ii);
    fprintf(fid, '\n Line Loop(O%d) = { ', Ostart+ii);
    for jj = 1:3
        
        if(Patch(ii,jj)<0)
            fprintf(fid, '-L%d, ', Lstart+abs(Patch(ii,jj)));
        else
            fprintf(fid, 'L%d, ', Lstart+abs(Patch(ii,jj)));
        end
        
    end
        
    if(Patch(ii,4)<0)
        fprintf(fid, '-L%d };', Lstart+abs(Patch(ii,4)));
    else
        fprintf(fid, 'L%d };', Lstart+abs(Patch(ii,4)));
    end
    
end


fprintf(fid, '\n\n  ');
% flag physical lines
for ii = 1:size(Ports,1)

    fprintf(fid, '\n\n Port%d = %d;',Nstart+2*ii-1, Nstart+2*ii-1);
    fprintf(fid, '\n Physical Line(Port%d) = {L%d};',Nstart+2*ii-1, Lstart+Ports(ii,1));
    fprintf(fid, '\n Port%d = %d;',Nstart+2*ii, Nstart+2*ii);
    fprintf(fid, '\n Physical Line(Port%d) = {L%d};',Nstart+2*ii, Lstart+Ports(ii,2));
    
end


fprintf(fid, '\n\n  ');
% print the surfaces 
for ii = 1:size(Patch,1)
    
    fprintf(fid, '\n S%d = %d;', Sstart+ii, Sstart+ii);
    fprintf(fid, '\n Ruled Surface(S%d) = {O%d};',Sstart+ii,Ostart+ii);
    
%     fprintf(fid, '\n Patch%d = %d;', Tstart+ii, Tstart+ii);
%     fprintf(fid, '\n Physical Surface(Patch%d) = {S%d};', Tstart+ii, Sstart+ii);
    
end


fprintf(fid, '\n\n  ');
% Physical Coil 
fprintf(fid, '\n Coil%d = %d;', Cnum, Nstart);
fprintf(fid, '\n Physical Surface(Coil%d) = {', Cnum);
for ii = 1:size(Patch,1)-1
    
    fprintf(fid, ' S%d,', Sstart+ii);
    
end
fprintf(fid, ' S%d };', Sstart+size(Patch,1));

fprintf(fid, '\n\n // --------------------------------------------------------------- '); 
fprintf(fid, '\n\n'); 

% -------------------------------------------------------------------------
% COIL DONE

fclose(fid);


% add Shield if exist
if (~isempty(ShieldRad)) && (~isempty(ShieldLen))
    [ok] = Shield_SCOIL_Gen(geofilename,ShieldRad,Xpos,Ypos,Zpos,ShieldLen);
end

% -------------------------------------------------------------------------
% call gmsh
% -------------------------------------------------------------------------
if ispc
    command = sprintf('.\\src_generate\\src_scoil\\gmsh %s', geofilename);
else
    if ismac
        command = sprintf('./src_generate/src_scoil/gmsh.app/Contents/MacOS/gmsh %s', geofilename);
    end
end
[ok, result] = system(command); 

