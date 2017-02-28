%%    MARIE installation script
% _________________________________________________________________________
%
%
%   Def.:   mex all functions needed in this project
%           First, choose your compiler with mex -setup
%
% _________________________________________________________________________
%
%     This program is part of the MARIE suite
%     MARIE - Magnetic Resonance Integral Equation suite
%     Copyright (C)2014 Jorge Fernandez Villena / Athanasios G. Polimeridis 
%     Computational Prototyping Group
%     Research Laboratory of Electronics
%     Massachusetts Institute of Technology
%     contact: jvillena@mit.edu / thanos_p@mit.edu
% 
%     MARIE is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     MARIE is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% _________________________________________________________________________
%%
[root,~,~] = fileparts(which('marie_install.m'));
cd(root);
clear root;

%%
try
	%% Cubatures
	cd(fullfile('.','src_vie','src_operators','mexcubatures'));
	mexdirect_ws_const_build;
	cd(current);

	%% coupling
	cd(fullfile('.','src_sie','src_assembly','src_mexdirect_ws_rwg'));
	mexdirect_ws_rwg_build;
	cd(current);
	cd(fullfile('.','src_sie','src_coupling','src_quadrature','coupling_mex'));
	mexCoupling_build;
	cd(current);
	cd(fullfile('.','src_sie','src_coupling','src_quadraturemex','coupling_qmex'));
	mexCoupling_build;
	cd(current);
	cd(fullfile('.','src_sie','src_coupling','src_quadraturemex','coupling_qmex'));
	mexCouplingACA_build;
	cd(current);

	%% data parsing
	cd(fullfile('.','src_geometry','src_scat','src_data'));
	mexDataBuild;
	cd(current);
catch me
	cd(current);
	error('Unable to compile MARIE\nMake sure to properly configure mex and openmp libraries.\nIdentifier: %s\nMessage: %s\n', me.identifier, me.message);
end
clear current;