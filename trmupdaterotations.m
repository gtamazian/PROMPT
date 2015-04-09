function result = trmupdaterotations(trmodel)
%TRMUPDATEROTATIONS Update rotation matrices by the optimal superposition
%   TRMUPDATEROTATIONS(trmodel) updates rotation matrices assigned to each
%   configuration of the specified transformation trmodel using the Kabsch
%   algorithm.
%
%   See also trmcreate
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

prevCoords = trmodel.StartCoords;
nModels = length(trmodel.U);

for j = 2:nModels
    currCoords = restorecoords(trmodel.r(:,j), trmodel.alpha(:,j), ...
        trmodel.psi(:,j));
    [~, prevCoords, transform] = procrustes(prevCoords, currCoords, ...
        'scaling', false, 'reflection', false);
    trmodel.U{j} = transform.T;
end

result = trmodel;

end

