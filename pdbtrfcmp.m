function result = pdbtrfcmp(PDBStruct1, PDBStruct2)
%PDBTRFCMP Compare transformations from PDB structures
%   PDBTRFCMP (PDBStruct1, PDBStruct2) compares transformations specified 
%   by PDB structures PDBStruct1 and PDBStruct2. The function returns a 
%   structure with the following fields:
%
%       * norm - norms of Cartesian coordinate difference matrices;
%       * rmsd - RMSDs after the Kabsch superposition;
%       * r - norms of bond length difference matrices;
%       * alpha - norms of planar angle circular difference matrices;
%       * psi - norms of torsion angle circular difference matrices.
%
%   Each field is a vector which elements correspond to a pair of models
%   from the specified structures.
%
%   See also createmodel pdbextractcoords
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

coords1 = pdbextractcoords(PDBStruct1);
coords2 = pdbextractcoords(PDBStruct2);

[r1, alpha1, psi1] = createmodel(PDBStruct1);
[r2, alpha2, psi2] = createmodel(PDBStruct2);

nModels = length(coords1); % the number of configurations

% get differences and RMSDs between configuration atom coordinates
normDiff = zeros(nModels, 1);
rmsd = zeros(nModels, 1);
rDiff = zeros(nModels, 1);
alphaDiff = zeros(nModels, 1);
psiDiff = zeros(nModels, 1);
for iModel = 1:nModels
    normDiff(iModel) = norm(coords1{iModel} - coords2{iModel});
    % superpose coords2{iModel} to coords1{iModel}
    [~, superposedCoords] = procrustes(coords1{iModel}, ...
        coords2{iModel}, 'scaling', false, 'reflection', false);
    rmsd(iModel) = mean(sqrt(sum((coords1{iModel} - ...
        superposedCoords).^2, 2)));
    rDiff(iModel) = norm(r1 - r2);
    alphaDiff(iModel) = norm(circdist(alpha1, alpha2));
    psiDiff(iModel) = norm(circdist(psi1, psi2));
end

result = struct('norm', normDiff, 'rmsd', rmsd, 'r', rDiff, ...
    'alpha', alphaDiff, 'psi', psiDiff);

end

