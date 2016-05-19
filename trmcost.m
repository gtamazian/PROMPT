function [cost, X, confX] = trmcost(trmodel)
%TRMCOST Calculate transformation cost
%   [cost, X, confX] = TRMCOST(trmodel) returns the cost of the 
%   transformation represented by the specified model. The cost is 
%   calculated as a function of the weighted sum of squared distances 
%   between atoms of the adjacent model configurations. Also the function 
%   returns restored coordinates of transformation atoms as an array of 
%   matrices X and restored coordinates of the unaligned configuration
%   atoms as an array of matrices confX.
%
%   See also trmobjfunc pdbtrfcost
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014-2016.
% mail (at) gtamazian (dot) com

p = 2;

[coords, confX] = trmrestorecoords(trmodel);

% calculate distances
distances = sqrt(sum((coords(:,:,2:end) - coords(:,:,1:end-1)).^2, 2));

cost = sum(trmodel.m .* sum(distances.^p, 3));
X = coords;

end

