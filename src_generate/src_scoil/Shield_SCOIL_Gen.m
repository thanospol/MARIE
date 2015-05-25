function [ok] = Shield_SCOIL_Gen(geofilename,Rad,Xpos,Ypos,Zpos,Len)
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
%       Rad - radius: the distance of the coil array to the center
%       Xpos,Ypos,Zpos - position of the center of coil array
%       Len - external length in z direction of each coil
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
if (Rad <= 0) 
    Rad = 0.2;
end

if (Len <= 0) 
    Len = 0.2;
end
  

% generate the name for the geometry file
fid = fopen(geofilename, 'a');
fprintf(fid, '\n\n\n // --------------------------------------------------------------- '); 
fprintf(fid, '\n // --------------------------------------------------------------- ');
fprintf(fid, '\n // CYLINDRIC SHIELD'); 
fprintf(fid, '\n // geometry created using MARIE coil generator');  
fprintf(fid, '\n // --------------------------------------------------------------- ');
fprintf(fid, '\n // --------------------------------------------------------------- '); 
Snum = 1000;

Arot = 0;
Res = 0.04;

pcount = 0;
lcount = 0;

Patch = zeros(4,4);

% -------------------------------------------------------------------------
% Top part of the coil

x = 0; y = 0; z = +Len/2; % center for top curve
Points = [x, y, z];
pcount = pcount + 1;
x = 0; y = 0; z = -Len/2; % center for top external bottom curve
Points = [Points; x, y, z];
pcount = pcount + 1;


% points ini --------------------------------------------------
Phi = 0;
x = Rad*sin(Phi); y = Rad*cos(Phi); z = +Len/2;
Points = [Points; x, y, z];
pcount = pcount + 1;
idxp1 = pcount;
x = Rad*sin(Phi); y = Rad*cos(Phi); z = -Len/2;
Points = [Points; x, y, z];
pcount = pcount + 1;
idxp2 = pcount;

Lines = [idxp1, idxp2, -100];
lcount = lcount + 1;
idxl1 = lcount;

for ii = 1:3

    % points end --------------------------------------------------
    Phi = ii*pi/2;
    x = Rad*sin(Phi); y = Rad*cos(Phi); z = -Len/2;
    Points = [Points; x, y, z];
    pcount = pcount + 1;
    idxp3 = pcount;
    x = Rad*sin(Phi); y = Rad*cos(Phi); z = +Len/2;
    Points = [Points; x, y, z];
    pcount = pcount + 1;
    idxp4 = pcount;
    
    
    % lines --------------------------------------------------
    Lines = [Lines; idxp2, idxp3, 2];
    lcount = lcount + 1;
    idxl2 = lcount;
    Lines = [Lines; idxp3, idxp4, -100];
    lcount = lcount + 1;
    idxl3 = lcount;
    Lines = [Lines; idxp4, idxp1, 1];
    lcount = lcount + 1;
    idxl4 = lcount;
    
    
    % surface --------------------------------------------------
    Patch(ii,:) = [idxl1 idxl2 idxl3 idxl4];


    % next patch
    idxp1 = idxp4;
    idxp2 = idxp3;
    idxl1 = -idxl3;
    
end


idxp3 = 4;
idxp4 = 3;
idxl3 = -1;

% lines --------------------------------------------------
Lines = [Lines; idxp2, idxp3, 2];
lcount = lcount + 1;
idxl2 = lcount;
Lines = [Lines; idxp4, idxp1, 1];
lcount = lcount + 1;
idxl4 = lcount;


% surface --------------------------------------------------
Patch(4,:) = [idxl1 idxl2 idxl3 idxl4];
    


% Print all the stuff in the file



% fprintf(fid, '\n Mesh.Format      = 1; // 1=.msh format by default');
% fprintf(fid, '\n Mesh.Algorithm   = 2; // 2D mesh algorithm (1=MeshAdapt, 5=Delaunay, 6=Frontal) Default value: 1 ');


fprintf(fid, '\n\n\n // --------------------------------------------------------------- '); 
fprintf(fid, '\n // Shield %d', Snum);  
fprintf(fid, '\n //        Radius       %g m', Rad);
fprintf(fid, '\n //        Translation Z  %g m', Zpos);
fprintf(fid, '\n //        Translation Y  %g m', Ypos);
fprintf(fid, '\n //        Translation X  %g m', Xpos);

fprintf(fid, '\n\n // Resolution');
fprintf(fid, '\n RES = %f;', Res);

Nstart = Snum*1000;
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
    fprintf(fid, '\n Point(P%d) = {%2.16g, %2.16g, %2.16g, RES};',Pstart+ii, x, y, z);
    
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
% print the surfaces 
for ii = 1:size(Patch,1)
    
    fprintf(fid, '\n S%d = %d;', Sstart+ii, Sstart+ii);
    fprintf(fid, '\n Ruled Surface(S%d) = {O%d};',Sstart+ii,Ostart+ii);
    
%     fprintf(fid, '\n Patch%d = %d;', Tstart+ii, Tstart+ii);
%     fprintf(fid, '\n Physical Surface(Patch%d) = {S%d};', Tstart+ii, Sstart+ii);
    
end



fprintf(fid, '\n\n  ');
% Physical Coil 
fprintf(fid, '\n Shield%d = %d;', Snum, Nstart);
fprintf(fid, '\n Physical Surface(Shield%d) = {', Snum);
for ii = 1:size(Patch,1)-1
    
    fprintf(fid, ' S%d,', Sstart+ii);
    
end
fprintf(fid, ' S%d };', Sstart+size(Patch,1));

fprintf(fid, '\n\n // --------------------------------------------------------------- '); 
fprintf(fid, '\n\n'); 

% -------------------------------------------------------------------------
% COIL DONE

fclose(fid);

% -------------------------------------------------------------------------
% call gmsh
% -------------------------------------------------------------------------
command = sprintf('.\\gmsh %s', geofilename);
[ok, result] = system(command); 

