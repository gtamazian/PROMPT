function initialPoint = trminitialpoint(trmodel, planarIndices, ...
    torsionIndices)
%TRMINITIALPOINT Get an initial point for the optimization procedure.
%   trminitialpoint(trmodel, planarIndices, torsionIndices) returns an
%   initial point containing values of the planar and torsion angles which
%   indices are specified in the planarIndices and torsionIndices
%   variables.
%
%   See also trmobjfunc
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2016.
% gaik (dot) tamazian (at) gmail (dot) com

initialPointP = trmodel.alpha(planarIndices, 2:end-1);
initialPointT = trmodel.psi(torsionIndices, 2:end-1);
initialPoint = reduceangles([initialPointP(:); initialPointT(:)]);

end

