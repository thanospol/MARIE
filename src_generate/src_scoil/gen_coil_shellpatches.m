function [ok] = gen_coil_shellpatches(geofilename,Cnum,Tports,Rports,Bports,Lports,Rad,Xtrans,Ytrans,Ztrans,Arot,Len,Asp,Wid,Gap)
%%    Surface Conformal Coil Generator
% _________________________________________________________________________
%
%
%   Def.: Fucntion 
%
% _________________________________________________________________________
%
%
%% INPUT
%
%
%% OUTPUT
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
 

fid = fopen(geofilename, 'a');
ok = 1;

Res = 3*Wid; % resolution

len = Len - 2*Wid; % internal length
wsp = asin(Wid/Rad); % width span in radians
asp = Asp - 2*wsp; % internal angle span in rads
gsp = asin(Gap/Rad);% gap span in rads

% ports in each side (Top, Right, Bottom, Left)
Tports = Tports*pi/180;
Bports = Bports*pi/180;

Rescorner = Res;
Respatch = Res;
Resport = Res;


Tports = sort(Tports, 'ascend');
Rports = sort(Rports, 'descend');
Bports = sort(Bports, 'descend');
Lports = sort(Lports, 'ascend');

% we will have N structures:
% point array
% line structure: one for each line or circle, 

pcount = 0;
lcount = 0;

% -------------------------------------------------------------------------
% Top part of the coil

x = 0; y = 0; z = +Len/2; % center for top curve
Points = [x, y, z, Res];
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


% -------------------------------------------------------------------------
% Top Left corner patch 

% points --------------------------------------------------
x = Rad*sin(-Asp/2); y = Rad*cos(-Asp/2); z = +len/2;
Points = [Points; x, y, z, Respatch];
pcount = pcount + 1;
idxp1 = pcount;
x = Rad*sin(-Asp/2); y = Rad*cos(-Asp/2); z = +Len/2;
Points = [Points; x, y, z, Rescorner];
pcount = pcount + 1;
idxp2 = pcount;
x = Rad*sin(-asp/2); y = Rad*cos(-asp/2); z = +Len/2;
Points = [Points; x, y, z, Respatch];
pcount = pcount + 1;
idxp3 = pcount;
x = Rad*sin(-asp/2); y = Rad*cos(-asp/2); z = +len/2;
Points = [Points; x, y, z, Rescorner];
pcount = pcount + 1;
idxp4 = pcount;

% lines --------------------------------------------------
Lines = [idxp1, idxp2, -100];
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
Patch = [idxl1 idxl2 idxl3 idxl4];


% -------------------------------------------------------------------------
% Top Right corner patch 

% points --------------------------------------------------
x = Rad*sin(asp/2); y = Rad*cos(asp/2); z = +len/2;
Points = [Points; x, y, z, Rescorner];
pcount = pcount + 1;
idxp1 = pcount;
x = Rad*sin(asp/2); y = Rad*cos(asp/2); z = +Len/2;
Points = [Points; x, y, z, Respatch];
pcount = pcount + 1;
idxp2 = pcount;
x = Rad*sin(Asp/2); y = Rad*cos(Asp/2); z = +Len/2;
Points = [Points; x, y, z, Rescorner];
pcount = pcount + 1;
idxp3 = pcount;
x = Rad*sin(Asp/2); y = Rad*cos(Asp/2); z = +len/2;
Points = [Points; x, y, z, Respatch];
pcount = pcount + 1;
idxp4 = pcount;

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
% Bottom Right corner patch 

% points --------------------------------------------------
x = Rad*sin(asp/2); y = Rad*cos(asp/2); z = -Len/2;
Points = [Points; x, y, z, Respatch];
pcount = pcount + 1;
idxp1 = pcount;
x = Rad*sin(asp/2); y = Rad*cos(asp/2); z = -len/2;
Points = [Points; x, y, z, Rescorner];
pcount = pcount + 1;
idxp2 = pcount;
x = Rad*sin(Asp/2); y = Rad*cos(Asp/2); z = -len/2;
Points = [Points; x, y, z, Respatch];
pcount = pcount + 1;
idxp3 = pcount;
x = Rad*sin(Asp/2); y = Rad*cos(Asp/2); z = -Len/2;
Points = [Points; x, y, z, Rescorner];
pcount = pcount + 1;
idxp4 = pcount;

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
% Bottom Left corner patch 

% points --------------------------------------------------
x = Rad*sin(-Asp/2); y = Rad*cos(-Asp/2); z = -Len/2;
Points = [Points; x, y, z, Rescorner];
pcount = pcount + 1;
idxp1 = pcount;
x = Rad*sin(-Asp/2); y = Rad*cos(-Asp/2); z = -len/2;
Points = [Points; x, y, z, Respatch];
pcount = pcount + 1;
idxp2 = pcount;
x = Rad*sin(-asp/2); y = Rad*cos(-asp/2); z = -len/2;
Points = [Points; x, y, z, Rescorner];
pcount = pcount + 1;
idxp3 = pcount;
x = Rad*sin(-asp/2); y = Rad*cos(-asp/2); z = -Len/2;
Points = [Points; x, y, z, Respatch];
pcount = pcount + 1;
idxp4 = pcount;

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
% Now start with internal patches

Ntp = length(Tports);
Nrp = length(Rports);
Nbp = length(Bports);
Nlp = length(Lports);
totalports = Ntp + Nrp + Nbp + Nlp;

Ports = zeros(totalports,2);
portcount = 0;
external = [1;2;6;7;11;12;16;13]; % external lines

% -------------------------------------------------------------------------
% Top patches

% start in the left top patch
% points  and lines defined
idxp1 = 8;
idxp2 = 7;
idxl1 = -3;

for ii = 1:Ntp
    
    % left points of the port ------------------------------------
    phi = Tports(ii) - gsp/2;
    x = Rad*sin(phi); y = Rad*cos(phi); z = +Len/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp3 = pcount;
    phi = Tports(ii) - gsp/2;
    x = Rad*sin(phi); y = Rad*cos(phi); z = +len/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp4 = pcount;
    
    % lines --------------------------------------------------
    Lines = [Lines; idxp2, idxp3, 1];
    lcount = lcount + 1;
    idxl2 = lcount;
    Lines = [Lines; idxp3, idxp4, -100];
    lcount = lcount + 1;
    idxl3 = lcount;
    Lines = [Lines; idxp4, idxp1, 2];
    lcount = lcount + 1;
    idxl4 = lcount;
    
    % port line ------------------------------------------------
    portcount = portcount + 1;
    Ports(portcount,1) = idxl3;
    external = [external; idxl2; idxl4];
    
    % surface --------------------------------------------------
    Patch = [Patch; idxl1 idxl2 idxl3 idxl4];
    
    % assign for the next patch
    
    % right points of the port ------------------------------------
    phi = Tports(ii) + gsp/2;
    x = Rad*sin(phi); y = Rad*cos(phi); z = +len/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp1 = pcount;
    phi = Tports(ii) + gsp/2;
    x = Rad*sin(phi); y = Rad*cos(phi); z = +Len/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp2 = pcount;
    
    % port line ------------------------------------------------
    Lines = [Lines; idxp1, idxp2, -100];
    lcount = lcount + 1;
    idxl1 = lcount;
    Ports(portcount,2) = idxl1;
    
end

% end in the right top patch
% points  and lines defined
idxp3 = 10;
idxp4 = 9;
idxl3 = -5;

Lines = [Lines; idxp2, idxp3, 1];
lcount = lcount + 1;
idxl2 = lcount;

Lines = [Lines; idxp4, idxp1, 2];
lcount = lcount + 1;
idxl4 = lcount;

% surface --------------------------------------------------
Patch = [Patch; idxl1, idxl2, idxl3, idxl4];
    


% -------------------------------------------------------------------------
% Right patches

% start in the top rigth patch
% points  and lines defined
idxp1 = 9;
idxp2 = 12;
idxl1 = -8;

for ii = 1:Nrp
    
    % top points of the port ------------------------------------
    x = Rad*sin(Asp/2); y = Rad*cos(Asp/2); z = Rports(ii) + Gap/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp3 = pcount;
    x = Rad*sin(asp/2); y = Rad*cos(asp/2); z = Rports(ii) + Gap/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp4 = pcount;
    
    % center point for curve line
    x = 0; y = 0; z = Rports(ii) + Gap/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    cp = pcount;
    
    % lines --------------------------------------------------
    Lines = [Lines; idxp2, idxp3, -100];
    lcount = lcount + 1;
    idxl2 = lcount;
    Lines = [Lines; idxp3, idxp4, cp];
    lcount = lcount + 1;
    idxl3 = lcount;
    Lines = [Lines; idxp4, idxp1, -100];
    lcount = lcount + 1;
    idxl4 = lcount;
    
    % port line ------------------------------------------------
    portcount = portcount + 1;
    Ports(portcount,1) = idxl3;
    external = [external; idxl2; idxl4];
    
    % surface --------------------------------------------------
    Patch = [Patch; idxl1 idxl2 idxl3 idxl4];
    
    % assign for the next patch
    
    % bottom points of the port ------------------------------------
    x = Rad*sin(asp/2); y = Rad*cos(asp/2); z = Rports(ii) - Gap/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp1 = pcount;
    x = Rad*sin(Asp/2); y = Rad*cos(Asp/2); z = Rports(ii) - Gap/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp2 = pcount;
    
    % center point for curve line
    x = 0; y = 0; z = Rports(ii) - Gap/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    cp = pcount;
    
   
    % port line ------------------------------------------------
    Lines = [Lines; idxp1, idxp2, cp];
    lcount = lcount + 1;
    idxl1 = lcount;
    Ports(portcount,2) = idxl1;
    
end

% end in the bottom right patch
% points  and lines defined
idxp3 = 15;
idxp4 = 14;
idxl3 = -10;

Lines = [Lines; idxp2, idxp3, -100];
lcount = lcount + 1;
idxl2 = lcount;

Lines = [Lines; idxp4, idxp1, -100];
lcount = lcount + 1;
idxl4 = lcount;

% surface --------------------------------------------------
Patch = [Patch; idxl1, idxl2, idxl3, idxl4];



    
% -------------------------------------------------------------------------
% Bottom patches

% start in the right bottom patch
% points  and lines defined
idxp1 = 14;
idxp2 = 13;
idxl1 = -9;

for ii = 1:Nbp
    
    % left points of the port ------------------------------------
    phi = Bports(ii) + gsp/2;
    x = Rad*sin(phi); y = Rad*cos(phi); z = -Len/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp3 = pcount;
    phi = Bports(ii) + gsp/2;
    x = Rad*sin(phi); y = Rad*cos(phi); z = -len/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp4 = pcount;
    
    % lines --------------------------------------------------
    Lines = [Lines; idxp2, idxp3, 4];
    lcount = lcount + 1;
    idxl2 = lcount;
    Lines = [Lines; idxp3, idxp4, -100];
    lcount = lcount + 1;
    idxl3 = lcount;
    Lines = [Lines; idxp4, idxp1, 3];
    lcount = lcount + 1;
    idxl4 = lcount;
    
    % port line ------------------------------------------------
    portcount = portcount + 1;
    Ports(portcount,1) = idxl3;
    external = [external; idxl2; idxl4];
    
    % surface --------------------------------------------------
    Patch = [Patch; idxl1 idxl2 idxl3 idxl4];
    
    % assign for the next patch
    
    % right points of the port ------------------------------------
    phi = Bports(ii) - gsp/2;
    x = Rad*sin(phi); y = Rad*cos(phi); z = -len/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp1 = pcount;
    phi = Bports(ii) - gsp/2;
    x = Rad*sin(phi); y = Rad*cos(phi); z = -Len/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp2 = pcount;
    
    % port line ------------------------------------------------
    Lines = [Lines; idxp1, idxp2, -100];
    lcount = lcount + 1;
    idxl1 = lcount;
    Ports(portcount,2) = idxl1;
    
end

% end in the left bottom patch
% points  and lines defined
idxp3 = 20;
idxp4 = 19;
idxl3 = -15;

Lines = [Lines; idxp2, idxp3, 4];
lcount = lcount + 1;
idxl2 = lcount;

Lines = [Lines; idxp4, idxp1, 3];
lcount = lcount + 1;
idxl4 = lcount;

% surface --------------------------------------------------
Patch = [Patch; idxl1, idxl2, idxl3, idxl4];




% -------------------------------------------------------------------------
% Left patches

% start in the bottom left patch
% points  and lines defined
idxp1 = 19;
idxp2 = 18;
idxl1 = -14;

for ii = 1:Nlp
    
    % top points of the port ------------------------------------
    x = Rad*sin(-Asp/2); y = Rad*cos(-Asp/2); z = Lports(ii) - Gap/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp3 = pcount;
    x = Rad*sin(-asp/2); y = Rad*cos(-asp/2); z = Lports(ii) - Gap/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp4 = pcount;
    
    % center point for curve line
    x = 0; y = 0; z = Lports(ii) - Gap/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    cp = pcount;
    
    % lines --------------------------------------------------
    Lines = [Lines; idxp2, idxp3, -100];
    lcount = lcount + 1;
    idxl2 = lcount;
    Lines = [Lines; idxp3, idxp4, cp];
    lcount = lcount + 1;
    idxl3 = lcount;
    Lines = [Lines; idxp4, idxp1, -100];
    lcount = lcount + 1;
    idxl4 = lcount;
    
    % port line ------------------------------------------------
    portcount = portcount + 1;
    Ports(portcount,1) = idxl3;
    external = [external; idxl2; idxl4];
    
    % surface --------------------------------------------------
    Patch = [Patch; idxl1 idxl2 idxl3 idxl4];
    
    % assign for the next patch
    
    % bottom points of the port ------------------------------------
    x = Rad*sin(-asp/2); y = Rad*cos(-asp/2); z = Lports(ii) + Gap/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp1 = pcount;
    x = Rad*sin(-Asp/2); y = Rad*cos(-Asp/2); z = Lports(ii) + Gap/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp2 = pcount;
    
    % center point for curve line
    x = 0; y = 0; z = Lports(ii) + Gap/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    cp = pcount;
    
   
    % port line ------------------------------------------------
    Lines = [Lines; idxp1, idxp2, cp];
    lcount = lcount + 1;
    idxl1 = lcount;
    Ports(portcount,2) = idxl1;
    
end

% end in the left top patch
% points  and lines defined
idxp3 = 5;
idxp4 = 8;
idxl3 = -4;

Lines = [Lines; idxp2, idxp3, -100];
lcount = lcount + 1;
idxl2 = lcount;

Lines = [Lines; idxp4, idxp1, -100];
lcount = lcount + 1;
idxl4 = lcount;

% surface --------------------------------------------------
Patch = [Patch; idxl1, idxl2, idxl3, idxl4];    



% Print all the stuff in the file


% fprintf(fid, '\n Mesh.Format      = 1; // 1=.msh format by default');
% fprintf(fid, '\n Mesh.Algorithm   = 2; // 2D mesh algorithm (1=MeshAdapt, 5=Delaunay, 6=Frontal) Default value: 1 ');


fprintf(fid, '\n\n\n // --------------------------------------------------------------- '); 
fprintf(fid, '\n // COIL %d', Cnum);  
fprintf(fid, '\n //        Radius       %g m', Rad);
fprintf(fid, '\n //        Rotation     %g rads', Arot);
fprintf(fid, '\n //        Translation Z  %g m', Ztrans);
fprintf(fid, '\n //        Translation Y  %g m', Ytrans);
fprintf(fid, '\n //        Translation X  %g m', Xtrans);

fprintf(fid, '\n\n // Resolution');
fprintf(fid, '\n RES = %f;', Res);

Nstart = Cnum*1000;
Pstart = Nstart + 100;
Lstart = Pstart + 400;
Ostart = Lstart + 200;
Sstart = Ostart + 100;
Tstart = Sstart + 100;

fprintf(fid, '\n\n  ');


% print the points
for ii = 1:pcount

    x = Points(ii,1)*cos(Arot) - Points(ii,2)*sin(Arot) + Xtrans;
    y = Points(ii,1)*sin(Arot) + Points(ii,2)*cos(Arot) + Ytrans;
    z = Points(ii,3) + Ztrans;    
    
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


    
