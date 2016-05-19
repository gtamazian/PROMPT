function result = trmchangeangles(trmodel, planarIndices, ...
    torsionIndices, angleValues )
%TRMCHANGEANGLES Update planar and torsion angle values in a model.
%   TRMCHANGEANGLES(trmodel, planarIndices, torsionIndices, angleValues)
%   updates values of the planar and torsion angles, which indices are
%   given in planarIndices and torsionIndices, respectively, with values
%   from angleValues.
%
%   See also trmdistantangleindices trmobjfunc
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2016.
% mail (at) gtamazian (dot) com

nConf = size(trmodel.r, 2);
nPlanar = length(planarIndices);
nTorsion = length(torsionIndices);

result = trmodel;

if ~isempty(planarIndices)
    result.alpha(planarIndices,2:end-1) = ...
        reshape(angleValues(1:nPlanar*(nConf - 2)), nPlanar, nConf - 2);
end

if ~isempty(torsionIndices)
    result.psi(torsionIndices,2:end-1) = ...
        reshape(angleValues(nPlanar*(nConf - 2)+1:end), nTorsion, ...
            nConf - 2);
end

end

