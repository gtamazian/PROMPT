function coords = pdbextractcoords(PDBStruct)
%PDBEXTRACTCOORDS Extract Cartesian atom coordinates from a PDB structure
%   PDBEXTRACTCOORDS(PDBStruct) extracts coordinates of the atoms that
%   constitute a protein from a PDB structure. The function returns a cell
%   array of matrices. Each matrix corresponds to a separate model from the
%   PDB structure.
%
%   See also restorecoords trmrestorecoords
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

nModels = length(PDBStruct.Model);
nAtoms = size(PDBStruct.Model(1).Atom, 2);
coords = zeros(nAtoms, 3, nModels);

for iModel = 1:nModels
    coords(:,:,iModel) = [[PDBStruct.Model(iModel).Atom.X]' ...
        [PDBStruct.Model(iModel).Atom.Y]' ...
        [PDBStruct.Model(iModel).Atom.Z]'];
end

end

