function result = trmupdaterotations(trmodel)
%TRMUPDATEROTATIONS Update rotation matrices by the optimal superposition
%   TRMUPDATEROTATIONS(trmodel) updates rotation matrices assigned to each
%   configuration of the specified transformation trmodel using the Kabsch
%   algorithm.
%
%   See also optimquat quat2rotmat
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

prevCoords = trmodel.StartCoords;
nModels = length(trmodel.U);

for j = 2:nModels
    currCoords = restorecoords(trmodel.r(:,j), trmodel.alpha(:,j), ...
        trmodel.psi(:,j));
    q = optimquat(prevCoords,currCoords);
    trmodel.U{j} = quat2rotmat(q);
    prevCoords = currCoords;
end

result = trmodel;

end

