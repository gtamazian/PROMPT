function U = trmsuperpos(trmodel)
% Calculate the rotation matrix for superposition of the second
% conformation to the first one.

nConf = size(trmodel.psi, 2) - 2;

U = cell(nConf+2, 1);
U{1} = eye(3);
prevCoords = trmodel.StartCoords;
for j = 2:nConf+2
    currCoords = restorecoords(trmodel.r(:,j), ...
        trmodel.alpha(:,j), trmodel.psi(:,j));
    [~,prevCoords,t] = procrustes(prevCoords, currCoords, ...
        'scaling', false, 'reflection', false);
    U{j} = t.T;
end

end
