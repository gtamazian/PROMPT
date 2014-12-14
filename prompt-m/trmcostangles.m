function cost = trmcostangles(trmodel, angleIndices, angles)
%TRMCOSTANGLES Transformation cost as a function of specified angles.
%   TRMCOSTANGLES(trmodel, angleIndices, angles) calculates cost of the
%   transformation defined by the model trmodel with the torsion angle 
%   values specified by the vector angle. Indices of the torsion angles 
%   are also specified in the vector angleIndices.
%
%   Example:
%
%       % Let t be a transformation model.
%       I = setdiff(1:size(t.psi, 1), 2:3:size(t.psi, 1));
%       trmcostangles(t, I, t.psi(I,2:end-1))
%
%   See also trmcost trmobjfunc pdbtrfcost
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

trmodel.psi(angleIndices,2:end-1) = ...
    reshape(angles, length(angleIndices), size(trmodel.psi, 2) - 2);
cost = trmcost(trmodel);

end

