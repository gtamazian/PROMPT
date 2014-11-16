function cost = trmcostangles(trmodel, indices, angles, p)
%TRMCOSTANGLES Transformation cost as a function of specified angles.
%   trmcostangles(trmodel, indices, angles) calculates cost of the trmodel
%	transformation with the torsion angle values specified by the
%   angle vector. Indices of the torsion angles are also specified.
%
%   Default values for the power p and cost type type are 2 and 'A',
%   respectively.
%
%   Example:
%       I = setdiff(1:size(t.psi, 1), 2:3:size(t.psi, 1));
%       trmcostangles(t, I, t.psi(I,2:end-1))
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

if nargin < 4
    p = 2;
end

trmodel.psi(indices,2:end-1) = ...
    reshape(angles, length(indices), size(trmodel.psi, 2) - 2);
cost = trmcost(trmodel, p);

end

