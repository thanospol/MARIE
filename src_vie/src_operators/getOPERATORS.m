function [op_out] = getOPERATORS(r,f,op,type,singular)
%% function that retrieves the operators N and K for JM-VIE
%
%% INPUT
% r         3-D Domain
% f         working frequency
%
% op        flag for choosing the operator
%                N: only operator N
%                K: only operator K
%
%% OPTIONAL INPUT
% type      the type for the chosen operator
%                empty: FFT of the Circulant (DEFAULT)
%                T: Toeplitz
%                C: Circulant
% singular  the method for the singular integrals in operator N
%                empty: DIRECTFN (DEFAULT)
%                DEMCEM
%
%% OUTPUT
% op_out    discrete operator
%
% -------------------------------------------------------------------------
%
%   A. G. Polimeridis -- thanos_p@mit.edu
%   J. Fernandez Villena -- jvillena@mit.edu
%   Computational Prototyping Group, RLE at MIT
%
% _________________________________________________________________________

% -------------------------------------------------------------------------
% Initialization of variables
% -------------------------------------------------------------------------

if(nargin < 4 || isempty(type))
   type = 'FFTofC';
end

if(nargin < 5 || isempty(singular))
   singular = 'DIRECTFN';
end

% -------------------------------------------------------------------------
% Domain & Frequency
% -------------------------------------------------------------------------
dx = r(2,1,1,1) - r(1,1,1,1); % voxel edge length
ko = 2*pi / 299792458 * f;    % wavenumber, ko = 2*pi/lambda = 2*pi/(co/f);
[Lr,Mr,Nr,~] = size(r);

% -------------------------------------------------------------------------
% Parallelization
% -------------------------------------------------------------------------
try
    if verLessThan('matlab', '8.3.0')
        isOpen = (matlabpool('size') > 0); %#ok<*DPOOL>
        if (~isOpen)
            mycluster = parcluster;
            matlabpool('local',mycluster.NumWorkers);
        end
    end
catch
    fprintf(1, '\n\n WARNING: matlabpool could not open, running on sequential');
    fprintf(1, '\n          computational performance may be compromised');
    fprintf(1, '\n          consider checking your MATLAB parallel capabilities');
    fprintf(1, '\n\n');
    pause(2);
end

fid = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          N OPERATOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(op ,'N')
    
    fprintf(fid, '\n ----------------------------------------------------------');
    fprintf(fid, '\n Generating the N Operator');
    fprintf(fid, '\n');
    
    % Singular integrals method
    if strcmp(singular ,'DIRECTFN')
        [opToeplitz] = assembly_nop(r,ko,dx,'DIRECTFN');
    elseif strcmp(singular ,'DEMCEM')
        [opToeplitz] = assembly_nop(r,ko,dx,'DEMCEM');
    end
    
    % Type
    if strcmp(type ,'FFTofC')
        opCirculant = circulant_nop(opToeplitz); 
        op_out = fft_operator(opCirculant); 
    elseif strcmp(type ,'T') 
        op_out = opToeplitz;           
    elseif strcmp(type ,'C')
        op_out = circulant_nop(opToeplitz);    
    end
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          K OPERATOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(op ,'K')
    
    fprintf(fid, '\n ----------------------------------------------------------');
    fprintf(fid, '\n Generating the K Operator');
    fprintf(fid, '\n');
    
    [opToeplitz] = assembly_kop(r,ko,dx);
    
    % Type
    if strcmp(type ,'FFTofC')
        opCirculant = circulant_kop(opToeplitz); 
        op_out = fft_operator(opCirculant); 
    elseif strcmp(type ,'T') 
        op_out = opToeplitz;           
    elseif strcmp(type ,'C')
        op_out = circulant_kop(opToeplitz);    
    end
    
end

infocir = whos('op_out');
fprintf(fid, '\n Frequency:           %3.3%f MHz', f/1e6);
fprintf(fid, '\n Domain Dimensions:   %dx%dx%d',Lr,Mr,Nr);
fprintf(fid, '\n FFT Dimensions:      %dx%dx%d',size(op_out,1),size(op_out,2), size(op_out,3));
fprintf(fid, '\n Memory space:        %.6f MB', infocir.bytes/(1024*1024));
fprintf(fid, '\n ----------------------------------------------------------\n ');

