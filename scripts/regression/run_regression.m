function [failed] = run_regression(UPDATE,FILENAME)
%%    Fucntion to run regression, and update or add new cases
% _________________________________________________________________________
%
%       Run a series of scripts and compare the solution with golden
%       if update  = 1 generates a new golden solution
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



% -------------------------------------------------------------------------
% Initialization of variables
% -------------------------------------------------------------------------


if(nargin < 1 || isempty(UPDATE))
   UPDATE = 0;
   fid = 1;
end

if(nargin < 2 || isempty(FILENAME))
    fid  = 1;
else
    fid = fopen(FILENAME, 'W');
    if fid  == 0
        fid  = 1;
        fprintf(fid, '\n\n\n Warning! Unable to open file for storing regression results: results to stdout\n\n');
    end
end


% -------------------------------------------------------------------------
% Add the corresponding path 
% -------------------------------------------------------------------------

% find the current folder
currentfolder = pwd;

% find the MARIE folder (two folders up)
if ispc
    idx = find(currentfolder == '\');
else
    idx = find(currentfolder == '/');
end
mariefolder = currentfolder(1:idx(end-1)-1);

% obtain the string with the recursive paths
p = genpath(mariefolder);

% add the path tree to the current path
addpath(p);



% -------------------------------------------------------------------------
% START with the REGRESSION
% -------------------------------------------------------------------------

REG_RESULTS = [];

% Add the regressions testcases here:


% -------------------------------------------------------------------------
%   S parameter of multi-port coil with homogeneous cube


[REG] = regression_svie_sparameters(UPDATE);
REG_RESULTS{1} = REG; 


% -------------------------------------------------------------------------
%   S parameter of multi-port coil with homogeneous cube


%[REG] = regression_svie_sparameters(UPDATE);
REG_RESULTS{2} = REG; 



% -------------------------------------------------------------------------
% REGRESSION DONE: REPORT!
% -------------------------------------------------------------------------

clc;
fprintf(fid, ' \n --------------------------------------------------------------');
fprintf(fid, ' \n --------------------------------------------------------------');
fprintf(fid, ' \n');
fprintf(fid, ' \n REGRESSION DONE... RESULTS:');
fprintf(fid, ' \n');

failed = 0;

for ii = 1:length(REG_RESULTS)
    
    
    REG = REG_RESULTS{ii};
    
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n --------------------------------------------------------------');
    fprintf(fid, ' \n TESTCASE: %s', REG.TESTCASE);
    fprintf(fid, ' \n');
    fprintf(fid, ' \n Elapsed time: %.2f [sec] (golden case %.2f [sec])', REG.TNEW,REG.TOLD);
    for jj = 1:length(REG.NAME)
        fprintf(fid, ' \n');
        fprintf(fid, ' \n RESULT: %s', REG.NAME{jj});
        fprintf(fid, ' \n ERROR:  %g', REG.ERROR(jj));
        if ( REG.ERROR(1) < 1e-12)
            fprintf(fid, ' \n         PASSED!');
        else
            fprintf(fid, ' \n         FAILED!');
            failed = failed + 1;
        end
    end

end


fprintf(fid, ' \n');
fprintf(fid, ' \n');
fprintf(fid, ' \n');
if failed
    fprintf(fid, ' \n REGRESSION FAILED! %d wrong  results', failed);
else
    fprintf(fid, ' \n REGRESSION PASSED! all cases passed');
end
fprintf(fid, ' \n');
fprintf(fid, ' \n --------------------------------------------------------------');
fprintf(fid, ' \n --------------------------------------------------------------');




