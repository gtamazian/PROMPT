function U = trmsuperpos(trmodel)
% Calculate the rotation matrix for superposition of the second
% conformation to the first one.

nConf = size(trmodel.psi, 2) - 2;
nAtoms = size(trmodel.r, 1) + 1;

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

U = cell(nConf+2, 1);
U{1} = eye(3);
prevCoords = trmodel.StartCoords(aligningAtoms,:);
for j = 2:nConf+2
    coords = restorecoords(trmodel.r(:,j), ...
        trmodel.alpha(:,j), trmodel.psi(:,j));
    currCoords = coords(aligningAtoms,:);
    [~,prevCoords,t] = procrustes(prevCoords, currCoords, ...
        'scaling', false, 'reflection', false);
    U{j} = t.T;
end

end
