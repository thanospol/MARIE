function [ZbcL,ZbcR] = Assembly_SCOUP_Q_BlockFactor(Scoord,index,etod,node,elem,freq,LEVEL_DVrule,tol,blockSize) 
%%    Quadrature coupling for the SIE+VIE solver
% _________________________________________________________________________
%
%   Fucntion to generate the Quadrature Coupling SIE to VIE
%   Goes triangle by triangle, and calls the mexCoupling function for each
%   PARFOR version
%   Generates a factorized coupling
%
% _________________________________________________________________________
%
%% Input
%       Scoord - coordinates of the observation points (No x 3)
%       node - coordinates of the nodes 
%       elem - 3 indexes of the nodes defining an element
%       etod - etod
%       index - mapping of the internal edge number to dof number
%       freq - frequency
%       LEVEL_DVrule - level for the DUNAVANT rule
%       tol - tolerance for the block SVD truncation
%
%
%% Output
%       ZbcL - 3xNo x Rank
%       ZbcR - Rank x Nd
%       Zbc ~ ZbcL*ZbcR in matrix form
%       Zbc - Tensor (No x 3 x Nd) with the contribution of each edge of the element
%              Zbc(1000,3,200) is z contribution of 200-th edge to 1000-th element
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
%            Define EM constants
% -------------------------------------------------------------------------

mu = 4*pi*1e-7;
co = 299792458;
eo = 1/co^2/mu;
%
omega = 2 * pi * freq;
lambda  = co/freq;
ko = 2*pi/lambda;

% Free-space impedance
% eta = omega*mu/ko; %3.767303134617706e+002; 

% -------------------------------------------------------------------------
% Define variables and allocate space
% -------------------------------------------------------------------------

NO = size(Scoord,1); % number of observation points
NE = size(elem,2); % number of elements
NC = max(index); % number of dofs

RO = Scoord.';

R1 = node(:,elem(1,:)); % 3xNe with coordinates of the first node of all elements
R2 = node(:,elem(2,:)); % 3xNe with coordinates of the first node of all elements
R3 = node(:,elem(3,:)); % 3xNe with coordinates of the first node of all elements

% get the sign and index for the contribution
ABSNUM = abs(etod(:,:)); % internal index of the edge
MULT = etod(:,:)./ABSNUM; % +1 or -1
IDX = index(ABSNUM);


% -------------------------------------------------------------------------
% 1D cubature's number of points
% -------------------------------------------------------------------------

% LEVEL_DVrule = 5; % usually more than enough
[ Np_2D, Z1, Z2, Z3, wp ] = dunavant_rule ( LEVEL_DVrule );

% [Np_2D,wp,Z1,Z2,Z3] = Gauss_2Dt(1);
% -------------------------------------------------------------------------
% loop on the elements and fill the matrix
% -------------------------------------------------------------------------

ZbcL = [];
ZbcR = [];
svTresh = [];

iniElem = 0;
blockSize = min(blockSize,NE-iniElem);

while (blockSize > 0)

%     Zblock = zeros(3*NO,blockSize*3);
%     for ii = 1:blockSize
%         % call the mex fun: for each element and all observation points
%         LOCALIDX = [3*ii-2; 3*ii-1; 3*ii];
%         Zblock = Zblock + mexQuadCoil2Scat(R1(:,ii+iniElem),R2(:,ii+iniElem),R3(:,ii+iniElem),1,RO(:),NO,LOCALIDX,MULT(:,ii+iniElem),3*blockSize,ko,Np_2D,Z1,Z2,Z3,wp);
%     end
    
    % corresponding blockSize elements
    % call the mex fun: for each element and all observation points
    LOCALIDX = 1:3*blockSize;
    LOCALELEM = iniElem+1:blockSize+iniElem;
    if (exist('ompQuadCoil2Scat', 'file') == 3)
        Zblock = ompQuadCoil2Scat(R1(:,LOCALELEM),R2(:,LOCALELEM),R3(:,LOCALELEM),blockSize,RO(:),NO,LOCALIDX,MULT(:,LOCALELEM),3*blockSize,ko,Np_2D,Z1,Z2,Z3,wp);
    else
        Zblock = mexQuadCoil2Scat(R1(:,LOCALELEM),R2(:,LOCALELEM),R3(:,LOCALELEM),blockSize,RO(:),NO,LOCALIDX,MULT(:,LOCALELEM),3*blockSize,ko,Np_2D,Z1,Z2,Z3,wp);
    end
    
    % backorth to previuos block
    if ~isempty(ZbcL)
        alpha = ZbcL'*Zblock;
        Zblock = Zblock - ZbcL*alpha;
    end
    
    % apply SVD on remaining block
    [Ublock, Sblock, Vblock] = svd(Zblock,'econ');
    clear Zblock;
    
    % define treshold for truncation
    if isempty(svTresh)
        svTresh = tol*Sblock(1,1);
    end
    % find the breaking point
    for kk = 1:length(Sblock)
        if (Sblock(kk,kk) < svTresh)
            break
        end
    end
    % truncate
    Ublock = Ublock(:,1:kk);
    Sblock = Sblock(1:kk,1:kk);
    Vblock = Vblock(:,1:kk);
    Vblock = Sblock*Vblock';
    clear Sblock;
    
    % store the new blocks in corresponding positions
    ZbcL = [ZbcL Ublock]; % append columns
    clear Ublock;
    
    % store the right factor in the corresponding positions
    % note that positions are dictated by index
    % and several elements may have the same index
    % so we have to take overlaps into account
    V = zeros(size(Vblock,1),NC);
    if ~isempty(ZbcR)
        for ii = 1:blockSize % in case there is overlapping (colIdx repeated)
            colIdx = IDX(1,ii+iniElem);
            if (colIdx) % we have existing 1st edge
                V(:,colIdx) = V(:,colIdx) + Vblock(:,3*ii-2);
                ZbcR(:,colIdx) = ZbcR(:,colIdx) + alpha(:,3*ii-2);
            end
            colIdx = IDX(2,ii+iniElem);
            if (colIdx) % we have existing 2nd edge
                V(:,colIdx) = V(:,colIdx) + Vblock(:,3*ii-1);
                ZbcR(:,colIdx) = ZbcR(:,colIdx) + alpha(:,3*ii-1);
            end
            colIdx = IDX(3,ii+iniElem);
            if (colIdx) % we have existing 3rd edge
                V(:,colIdx) = V(:,colIdx) + Vblock(:,3*ii);
                ZbcR(:,colIdx) = ZbcR(:,colIdx) + alpha(:,3*ii);
            end
        end
    else
        for ii = 1:blockSize % in case there is overlapping (colIdx repeated)
            colIdx = IDX(1,ii+iniElem);
            if (colIdx) % we have existing 1st edge
                V(:,colIdx) = V(:,colIdx) + Vblock(:,3*ii-2);
            end
            colIdx = IDX(2,ii+iniElem);
            if (colIdx) % we have existing 2nd edge
                V(:,colIdx) = V(:,colIdx) + Vblock(:,3*ii-1);
            end
            colIdx = IDX(3,ii+iniElem);
            if (colIdx) % we have existing 3rd edge
                V(:,colIdx) = V(:,colIdx) + Vblock(:,3*ii);
            end
        end
    end
    clear alpha; clear Vblock;
    ZbcR = [ ZbcR; V]; % append rows
    
    % restart counters and arrays
    iniElem = iniElem + blockSize;
    blockSize = min(blockSize,NE-iniElem);
    
end
  
% Zbc = Zbc/(1j*ko);
% scalefactor = (1j*eta)/(4*pi*ko);
% -------------------------------------------------------------------------
%             Final Z (with mult. constant) 
%
%      4*pi comes from omitting it in Green function!!
%       Z = discretization of  -e^scattered
% -------------------------------------------------------------------------
ce = 1i*omega*eo;
scalefactor = - 1 / ce / (4*pi);

ZbcL = scalefactor * ZbcL;
%  Zbc = (4*pi)/(1i*eta) * Zbc; 
    