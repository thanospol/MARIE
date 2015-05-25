function [R_faces] = cube_faces(dx,dy,dz)
% 6 faces of cube
% S_-x
R_mins_x(:,1) = [0,0 ,0];
R_mins_x(:,2) = [0,dy,0];
R_mins_x(:,3) = [0,dy,dz];
R_mins_x(:,4) = [0,0 ,dz];
% S_+x
R_plus_x(:,1) = [dx,0 ,0];
R_plus_x(:,2) = [dx,dy,0];
R_plus_x(:,3) = [dx,dy,dz];
R_plus_x(:,4) = [dx,0 ,dz];
% S_-y
R_mins_y(:,1) = [0, 0,0];
R_mins_y(:,2) = [dx,0,0];
R_mins_y(:,3) = [dx,0,dz];
R_mins_y(:,4) = [0 ,0,dz];
% S_+y
R_plus_y(:,1) = [0 ,dy,0];
R_plus_y(:,2) = [dx,dy,0];
R_plus_y(:,3) = [dx,dy,dz];
R_plus_y(:,4) = [0 ,dy,dz];
% S_-z
R_mins_z(:,1) = [0, 0 ,0];
R_mins_z(:,2) = [dx,0 ,0];
R_mins_z(:,3) = [dx,dy,0];
R_mins_z(:,4) = [0 ,dy,0];
% S_+z
R_plus_z(:,1) = [0 ,0 ,dz];
R_plus_z(:,2) = [dx,0 ,dz];
R_plus_z(:,3) = [dx,dy,dz];
R_plus_z(:,4) = [0 ,dy,dz];
% Output
R_faces(:,:,1) = R_mins_x;
R_faces(:,:,2) = R_plus_x;
R_faces(:,:,3) = R_mins_y;
R_faces(:,:,4) = R_plus_y;
R_faces(:,:,5) = R_mins_z;
R_faces(:,:,6) = R_plus_z;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%