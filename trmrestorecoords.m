function coords = trmrestorecoords(trmodel)
%TRMRESTORECOORDS Restore Cartesian coordinates of transformation atoms
%   TRMRESTORECOORDS(trmodel) returns a cell array of matrices containing
%   Cartesian coordinates of the atoms that constitute configurations of
%   the transformation trmodel.
%
%   See also restorecoords
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

nConf = size(trmodel.psi, 2);
nAtoms = size(trmodel.r, 1) + 1;
coords = cell(1, nConf);
coords{1} = trmodel.StartCoords;

if ~isfield(trmodel, 'proteinType')
    aligningAtoms = 1:nAtoms;
elseif strcmp(trmodel.proteinType, 'loop')
    aligningAtoms = [1:3, nAtoms-2:nAtoms];
elseif strcmp(trmodel.proteinType, 'n-tail')
    aligningAtoms = 1:3;
elseif strcmp(trmodel.proteinType, 'c-tail')
    aligningAtoms = nAtoms-2:nAtoms;
else
    aligningAtoms = 1:nAtoms;
end

firstConfTranslation = repmat(mean(coords{1}(aligningAtoms,:),1),nAtoms,1);

for i = 2:nConf
    coords{i} = restorecoords(trmodel.r(:,i), ...
        trmodel.alpha(:,i), trmodel.psi(:,i));
    coords{i} = coords{i}*trmodel.U{i};
    
    % apply the translation
    currTranslation = repmat(mean(coords{i}(aligningAtoms,:), 1),nAtoms,1);
    coords{i} = coords{i} - currTranslation + firstConfTranslation;       
end

end

