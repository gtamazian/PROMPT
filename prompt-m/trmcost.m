function [cost, X] = trmcost(trmodel, p)
%TRMCOST Calculates cost of the transformation.
%   trmcost(trmodel) returns the cost of the transformation represented by
%   the specified model. The cost is calculated as a function of the
%   weighted sum of distances between atoms of the adjacent models. The
%   parameter p is the power interatomic distances are raised to.
%
%   The default value for the power p is 2.
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

if nargin < 2
    p = 2;
end

coordscell = trmrestorecoords(trmodel);
nModels = length(coordscell);
nAtoms = size(coordscell{1}, 1);

% calculate distances
distances = zeros(nAtoms, nModels);
for i = 2:nModels
   distances(:,i) = sqrt(sum((coordscell{i} - coordscell{i-1}).^2, 2));
end

cost = sum(trmodel.m .* sum(distances.^p, 2));
X = coordscell;

end

