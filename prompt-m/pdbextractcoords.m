function coords = pdbextractcoords(PDBStruct)
%PDBEXTRACTCOORDS Extracts atom coordinates from a PDB structure.
%   pdbextractcoords(PDBStruct) extracts coordinates of the atoms that
%   represent a protein from a PDB structure. A cell array of matrices is
%   returned; each matrix corresponds to a separate model.
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

nModels = length(PDBStruct.Model);
coords = cell(1, nModels);

for i = 1:nModels
    coords{i} = [[PDBStruct.Model(i).Atom.X]' ...
        [PDBStruct.Model(i).Atom.Y]' [PDBStruct.Model(i).Atom.Z]'];
end

end

