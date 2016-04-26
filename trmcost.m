function [cost, X] = trmcost(trmodel)
%TRMCOST Calculate transformation cost
%   [cost, X] = TRMCOST(trmodel) returns the cost of the transformation 
%   represented by the specified model. The cost is calculated as a 
%   function of the weighted sum of squared distances between atoms of the 
%   adjacent model configurations. Also the function returns restored
%   coordinates of transformation atoms as a cell array of matrices X.
%
%   See also trmobjfunc pdbtrfcost
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

p = 2;

coords = trmrestorecoords(trmodel);
nModels = size(coords, 3);
nAtoms = size(coords(:,:,1), 1);

% calculate distances
distances = zeros(nAtoms, nModels);
for i = 2:nModels
   distances(:,i) = sqrt(sum((coords(:,:,i) - coords(:,:,i-1)).^2, 2));
end

cost = sum(trmodel.m .* sum(distances.^p, 2));
X = coords;

end

