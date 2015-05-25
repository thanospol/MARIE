function [figidx,fighandle] = plot_3Dvector(vec,r,xcut,ycut,zcut,figidx,scale,name)
%%   Plots a 3D vector for a given cut 
% _________________________________________________________________________
%
%       Plots the 3 maps of a given 3D vector in the same figure
%
% _________________________________________________________________________
%
%% INPUT
%   vec:    3D vector
%   r:      4D domain grid
%   xcut:   x coordinate for maps
%   ycut:   y coordinate for maps
%   zcut:   z coordinate for maps
%   figidx: figure count
%   scale:  array 2x3 with minimum and maximun for each map
%   name:   name of the vector for the plot
%
%
%% OUTPUT
%   figidx: new figure count
%   fighandle: handle to the figure
%
% -------------------------------------------------------------------------
%
%   J. Fernandez Villena -- jvillena@mit.edu
%   A.G. Polimeridis -- thanos_p@mit.edu
%   Computational Prototyping Group, RLE at MIT
%
% _________________________________________________________________________


% -------------------------------------------------------------------------
% Prepare data
% -------------------------------------------------------------------------

if (nargin < 6) || isempty(figidx)
    figidx = 0;
end
if (nargin < 7) || isempty(scale)
    scale = [];
else
    if length(scale)==1
        scale = scale*ones(3,2);
        scale(:,2) = 0;
    end
end
if (nargin < 8) || isempty(name)
    name = ' ';
end


[L,M,N,~] = size(r);

% find indexes for the cuts
[~, xidx] = min(abs(r(:,1,1,1) - xcut));
[~, yidx] = min(abs(r(1,:,1,2) - ycut));
[~, zidx] = min(abs(r(1,1,:,3) - zcut));


% -------------------------------------------------------------------------
% Plot figures
% -------------------------------------------------------------------------

figidx = figidx + 1;
figure(figidx);

% -------------------------------------------------------------------------
% Subplot 1

subplot(1,3,1);
if (~isempty(scale))
    Scaleplot = scale(1,1)*ones(M,N);
    Scaleplot(end,end) = scale(1,2);
    imagesc(rot90(Scaleplot));
    hold on;
end
imagesc(rot90(squeeze(vec(xidx,end:-1:1,:))));
set(gca, 'XTick', []);
set(gca, 'YTick', []);
title(sprintf('(sag.) %s', name), 'FontSize', 18);
colormap('hot');
colorbar;
axis image
set(gca, 'FontSize', 18);

% -------------------------------------------------------------------------
% Subplot 2

subplot(1,3,2);
if (~isempty(scale))
    Scaleplot = scale(2,1)*ones(L,N);
    Scaleplot(end,end) = scale(2,2);
    imagesc(rot90(Scaleplot));
    hold on;
end
imagesc(rot90(squeeze(vec(:,yidx,:))));
set(gca, 'XTick', []);
set(gca, 'YTick', []);
title(sprintf('(cor.) %s', name), 'FontSize', 18);
colormap('hot');
colorbar;
axis image
set(gca, 'FontSize', 18);


% -------------------------------------------------------------------------
% Subplot 3

subplot(1,3,3);
if (~isempty(scale))
    Scaleplot = scale(3,1)*ones(L,M);
    Scaleplot(end,end) = scale(3,2);
    imagesc(rot90(Scaleplot));
    hold on;
end
imagesc(rot90(squeeze(vec(:,:,zidx))));
set(gca, 'XTick', []);
set(gca, 'YTick', []);
title(sprintf('(ax.) %s', name), 'FontSize', 18);
colormap('hot');
colorbar;
axis image
set(gca, 'FontSize', 18);

fighandle = gcf;

