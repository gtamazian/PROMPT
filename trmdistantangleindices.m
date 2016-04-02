function angleIndices = trmdistantangleindices(trmodel, nAngles, angleType)
%TRMDISTANTANGLEINDICES Get indices of the most distant torsion angles
%   TRMDISTANTANGLEINDICES(trmodel, nAngles, angleType) returns indices of 
%   nAngles angles of the type specified by angleType which values in the 
%   first and last model configurations differ at most. The angleType
%   argument must be 'torsion' (the default value) or 'planar'.
%
%   See also trmobjfunc
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014-2016.
% gaik (dot) tamazian (at) gmail (dot) com

if nargin < 3
    angleType = 'torsion';
end

if strcmp(angleType, 'torsion')
    angleValues = trmodel.psi;
elseif strcmp(angleType, 'planar')
    angleValues = trmodel.alpha;
else
    error('PROMPT:trmplotanglediff:trmdistantangleindices', ...
        'Incorrect angle type specified.');
end

[~,angleIndices] = sort(abs(circdist(angleValues(:,1), ...
    angleValues(:,end))), 'descend');
angleIndices = angleIndices(1:nAngles);

end

