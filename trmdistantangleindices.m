function angleIndices = trmdistantangleindices(trmodel, nAngles)
%TRMDISTANTANGLEINDICES Get indices of the most distant torsion angles
%   TRMDISTANTANGLEINDICES(trmodel, nAngles) returns indices of nAngles
%   torsion angles which values in the first and last model configurations
%   differ at most.
%
%   See also trmobjfunc
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

[~,angleIndices] = sort(abs(circdist(trmodel.psi(:,1), ...
    trmodel.psi(:,end))), 'descend');
angleIndices = angleIndices(1:nAngles);

end

