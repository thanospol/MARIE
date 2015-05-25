function [ok] = gen_coil_planarpatches(geofilename,Cnum,Tports,Rports,Bports,Lports,Rad,Xtrans,Ytrans,Ztrans,Arot,Len,Ssp,Wid,Gap)
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
%   J. Fernandez Villena -- jvillena@mit.edu
%   A.G. Polimeridis -- thanos_p@mit.edu
%   Computational Prototyping Group, RLE at MIT
%
% _________________________________________________________________________

fid = fopen(geofilename, 'a');
ok = 1;

Res = 3*Wid;

len = Len - 2*Wid; % internal length
ssp = Ssp - 2*Wid; % internal side


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
% Top Left corner patch 

% points --------------------------------------------------
x = -Ssp/2; y = Rad; z = +len/2;
Points = [x, y, z, Respatch];
pcount = pcount + 1;
idxp1 = pcount;
x = -Ssp/2; y = Rad; z = +Len/2;
Points = [Points; x, y, z, Rescorner];
pcount = pcount + 1;
idxp2 = pcount;
x = -ssp/2; y = Rad; z = +Len/2;
Points = [Points; x, y, z, Respatch];
pcount = pcount + 1;
idxp3 = pcount;
x = -ssp/2; y = Rad; z = +len/2;
Points = [Points; x, y, z, Rescorner];
pcount = pcount + 1;
idxp4 = pcount;

% lines --------------------------------------------------
Lines = [idxp1, idxp2];
lcount = lcount + 1;
idxl1 = lcount;
Lines = [Lines; idxp2, idxp3];
lcount = lcount + 1;
idxl2 = lcount;
Lines = [Lines; idxp3, idxp4];
lcount = lcount + 1;
idxl3 = lcount;
Lines = [Lines; idxp4, idxp1];
lcount = lcount + 1;
idxl4 = lcount;

% surface --------------------------------------------------
Patch = [idxl1 idxl2 idxl3 idxl4];


% -------------------------------------------------------------------------
% Top Right corner patch 

% points --------------------------------------------------
x = ssp/2; y = Rad; z = +len/2;
Points = [Points; x, y, z, Rescorner];
pcount = pcount + 1;
idxp1 = pcount;
x = ssp/2; y = Rad; z = +Len/2;
Points = [Points; x, y, z, Respatch];
pcount = pcount + 1;
idxp2 = pcount;
x = Ssp/2; y = Rad; z = +Len/2;
Points = [Points; x, y, z, Rescorner];
pcount = pcount + 1;
idxp3 = pcount;
x = Ssp/2; y = Rad; z = +len/2;
Points = [Points; x, y, z, Respatch];
pcount = pcount + 1;
idxp4 = pcount;

% lines --------------------------------------------------
Lines = [Lines; idxp1, idxp2];
lcount = lcount + 1;
idxl1 = lcount;
Lines = [Lines; idxp2, idxp3];
lcount = lcount + 1;
idxl2 = lcount;
Lines = [Lines; idxp3, idxp4];
lcount = lcount + 1;
idxl3 = lcount;
Lines = [Lines; idxp4, idxp1];
lcount = lcount + 1;
idxl4 = lcount;

% surface --------------------------------------------------
Patch = [Patch; idxl1 idxl2 idxl3 idxl4];


% -------------------------------------------------------------------------
% Bottom Right corner patch 

% points --------------------------------------------------
x = ssp/2; y = Rad; z = -Len/2;
Points = [Points; x, y, z, Respatch];
pcount = pcount + 1;
idxp1 = pcount;
x = ssp/2; y = Rad; z = -len/2;
Points = [Points; x, y, z, Rescorner];
pcount = pcount + 1;
idxp2 = pcount;
x = Ssp/2; y = Rad; z = -len/2;
Points = [Points; x, y, z, Respatch];
pcount = pcount + 1;
idxp3 = pcount;
x = Ssp/2; y = Rad; z = -Len/2;
Points = [Points; x, y, z, Rescorner];
pcount = pcount + 1;
idxp4 = pcount;

% lines --------------------------------------------------
Lines = [Lines; idxp1, idxp2];
lcount = lcount + 1;
idxl1 = lcount;
Lines = [Lines; idxp2, idxp3];
lcount = lcount + 1;
idxl2 = lcount;
Lines = [Lines; idxp3, idxp4];
lcount = lcount + 1;
idxl3 = lcount;
Lines = [Lines; idxp4, idxp1];
lcount = lcount + 1;
idxl4 = lcount;

% surface --------------------------------------------------
Patch = [Patch; idxl1 idxl2 idxl3 idxl4];


% -------------------------------------------------------------------------
% Bottom Left corner patch 

% points --------------------------------------------------
x = -Ssp/2; y = Rad; z = -Len/2;
Points = [Points; x, y, z, Rescorner];
pcount = pcount + 1;
idxp1 = pcount;
x = -Ssp/2; y = Rad; z = -len/2;
Points = [Points; x, y, z, Respatch];
pcount = pcount + 1;
idxp2 = pcount;
x = -ssp/2; y = Rad; z = -len/2;
Points = [Points; x, y, z, Rescorner];
pcount = pcount + 1;
idxp3 = pcount;
x = -ssp/2; y = Rad; z = -Len/2;
Points = [Points; x, y, z, Respatch];
pcount = pcount + 1;
idxp4 = pcount;

% lines --------------------------------------------------
Lines = [Lines; idxp1, idxp2];
lcount = lcount + 1;
idxl1 = lcount;
Lines = [Lines; idxp2, idxp3];
lcount = lcount + 1;
idxl2 = lcount;
Lines = [Lines; idxp3, idxp4];
lcount = lcount + 1;
idxl3 = lcount;
Lines = [Lines; idxp4, idxp1];
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
idxp1 = 4;
idxp2 = 3;
idxl1 = -3;

for ii = 1:Ntp
    
    % left points of the port ------------------------------------
    x = Tports(ii) - Gap/2; y = Rad; z = +Len/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp3 = pcount;
    x = Tports(ii) - Gap/2; y = Rad; z = +len/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp4 = pcount;
    
    % lines --------------------------------------------------
    Lines = [Lines; idxp2, idxp3];
    lcount = lcount + 1;
    idxl2 = lcount;
    Lines = [Lines; idxp3, idxp4];
    lcount = lcount + 1;
    idxl3 = lcount;
    Lines = [Lines; idxp4, idxp1];
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
    x = Tports(ii) + Gap/2; y = Rad; z = +len/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp1 = pcount;
    x = Tports(ii) + Gap/2; y = Rad; z = +Len/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp2 = pcount;
    
    % port line ------------------------------------------------
    Lines = [Lines; idxp1, idxp2];
    lcount = lcount + 1;
    idxl1 = lcount;
    Ports(portcount,2) = idxl1;
    
end

% end in the right top patch
% points  and lines defined
idxp3 = 6;
idxp4 = 5;
idxl3 = -5;

Lines = [Lines; idxp2, idxp3];
lcount = lcount + 1;
idxl2 = lcount;

Lines = [Lines; idxp4, idxp1];
lcount = lcount + 1;
idxl4 = lcount;

% surface --------------------------------------------------
Patch = [Patch; idxl1, idxl2, idxl3, idxl4];
    


% -------------------------------------------------------------------------
% Right patches

% start in the top rigth patch
% points  and lines defined
idxp1 = 5;
idxp2 = 8;
idxl1 = -8;

for ii = 1:Nrp
    
    % top points of the port ------------------------------------
    x = Ssp/2; y = Rad; z = Rports(ii) + Gap/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp3 = pcount;
    x = ssp/2; y = Rad; z = Rports(ii) + Gap/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp4 = pcount;
    
   
    % lines --------------------------------------------------
    Lines = [Lines; idxp2, idxp3];
    lcount = lcount + 1;
    idxl2 = lcount;
    Lines = [Lines; idxp3, idxp4];
    lcount = lcount + 1;
    idxl3 = lcount;
    Lines = [Lines; idxp4, idxp1];
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
    x = ssp/2; y = Rad; z = Rports(ii) - Gap/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp1 = pcount;
    x = Ssp/2; y = Rad; z = Rports(ii) - Gap/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp2 = pcount;
    
     
   
    % port line ------------------------------------------------
    Lines = [Lines; idxp1, idxp2];
    lcount = lcount + 1;
    idxl1 = lcount;
    Ports(portcount,2) = idxl1;
    
end

% end in the bottom right patch
% points  and lines defined
idxp3 = 11;
idxp4 = 10;
idxl3 = -10;

Lines = [Lines; idxp2, idxp3];
lcount = lcount + 1;
idxl2 = lcount;

Lines = [Lines; idxp4, idxp1];
lcount = lcount + 1;
idxl4 = lcount;

% surface --------------------------------------------------
Patch = [Patch; idxl1, idxl2, idxl3, idxl4];



    
% -------------------------------------------------------------------------
% Bottom patches

% start in the right bottom patch
% points  and lines defined
idxp1 = 10;
idxp2 = 9;
idxl1 = -9;

for ii = 1:Nbp
    
    % left points of the port ------------------------------------
    x = Bports(ii) + Gap/2; y = Rad; z = -Len/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp3 = pcount;
    x = Bports(ii) + Gap/2; y = Rad; z = -len/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp4 = pcount;
    
    % lines --------------------------------------------------
    Lines = [Lines; idxp2, idxp3];
    lcount = lcount + 1;
    idxl2 = lcount;
    Lines = [Lines; idxp3, idxp4];
    lcount = lcount + 1;
    idxl3 = lcount;
    Lines = [Lines; idxp4, idxp1];
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
    x = Bports(ii) - Gap/2; y = Rad; z = -len/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp1 = pcount;
    x = Bports(ii) - Gap/2; y = Rad; z = -Len/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp2 = pcount;
    
    % port line ------------------------------------------------
    Lines = [Lines; idxp1, idxp2];
    lcount = lcount + 1;
    idxl1 = lcount;
    Ports(portcount,2) = idxl1;
    
end

% end in the left bottom patch
% points  and lines defined
idxp3 = 16;
idxp4 = 15;
idxl3 = -15;

Lines = [Lines; idxp2, idxp3];
lcount = lcount + 1;
idxl2 = lcount;

Lines = [Lines; idxp4, idxp1];
lcount = lcount + 1;
idxl4 = lcount;

% surface --------------------------------------------------
Patch = [Patch; idxl1, idxl2, idxl3, idxl4];




% -------------------------------------------------------------------------
% Left patches

% start in the bottom left patch
% points  and lines defined
idxp1 = 15;
idxp2 = 14;
idxl1 = -14;

for ii = 1:Nlp
    
    % top points of the port ------------------------------------
    x = -Ssp/2; y = Rad; z = Lports(ii) - Gap/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp3 = pcount;
    x = -ssp/2; y = Rad; z = Lports(ii) - Gap/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp4 = pcount;
    
    % lines --------------------------------------------------
    Lines = [Lines; idxp2, idxp3];
    lcount = lcount + 1;
    idxl2 = lcount;
    Lines = [Lines; idxp3, idxp4];
    lcount = lcount + 1;
    idxl3 = lcount;
    Lines = [Lines; idxp4, idxp1];
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
    x = -ssp/2; y = Rad; z = Lports(ii) + Gap/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp1 = pcount;
    x = -Ssp/2; y = Rad; z = Lports(ii) + Gap/2;
    Points = [Points; x, y, z, Resport];
    pcount = pcount + 1;
    idxp2 = pcount;
    
  
    % port line ------------------------------------------------
    Lines = [Lines; idxp1, idxp2];
    lcount = lcount + 1;
    idxl1 = lcount;
    Ports(portcount,2) = idxl1;
    
end

% end in the left top patch
% points  and lines defined
idxp3 = 1;
idxp4 = 4;
idxl3 = -4;

Lines = [Lines; idxp2, idxp3];
lcount = lcount + 1;
idxl2 = lcount;

Lines = [Lines; idxp4, idxp1];
lcount = lcount + 1;
idxl4 = lcount;

% surface --------------------------------------------------
Patch = [Patch; idxl1, idxl2, idxl3, idxl4];    



% Print all the stuff in the file


% fprintf(fid, '\n Mesh.Format      = 1; // 1=.msh format by default');
% fprintf(fid, '\n Mesh.Algorithm   = 2; // 2D mesh algorithm (1=MeshAdapt, 5=Delaunay, 6=Frontal) Default value: 1 ');


fprintf(fid, '\n\n\n // --------------------------------------------------------------- '); 
fprintf(fid, '\n // COIL %d', Cnum);  
fprintf(fid, '\n //        Radius         %g m', Rad);
fprintf(fid, '\n //        Rotation       %g rads', Arot);
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

    x = Points(ii,1)*cos(Arot) - Points(ii,2)*sin(Arot);
    y = Points(ii,1)*sin(Arot) + Points(ii,2)*cos(Arot);
    x = x + Xtrans;
    y = y + Ytrans;
    z = Points(ii,3) + Ztrans;
    
    fprintf(fid, '\n P%d = %d;',Pstart+ii, Pstart+ii);
    fprintf(fid, '\n Point(P%d) = {%2.16g, %2.16g, %2.16g, RES};',Pstart+ii, x, y, z);
    
end


fprintf(fid, '\n\n  ');

% print the lines
for ii = 1:lcount
         
        fprintf(fid, '\n L%d = %d;',Lstart+ii, Lstart+ii);
        fprintf(fid, '\n Line(L%d) = {P%d, P%d};',Lstart+ii, Pstart+Lines(ii,1), Pstart+Lines(ii,2));
   
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


    
