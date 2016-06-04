function outputStruct = pdbcartinterp(PDBStruct1, PDBStruct2, nConf)
%PDBCARTINTERP Interpolation of Cartesian atom coordinates
%   PDBCARTINTERP(PDBStruct1, PDBSTruct2, nConf) procudes a PDB structure
%   that contains nConf models which atom coordinates were obtained by
%   linear interpolation of their Cartesian coordinates.
%
%   See also pdbextractcoords trmcreate
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2016.
% mail (at) gtamazian (dot) com

coords1 = pdbextractcoords(PDBStruct1);
coords2 = pdbextractcoords(PDBStruct2);

[~, coords2] = procrustes(coords1, coords2, 'scaling', false, ...
    'reflection', false);

nAtoms = size(coords1, 1);

outputStruct = PDBStruct1;

interpMask1 = reshape(repelem(...
    repelem(linspace(0, 1, nConf), 3, 1), 1, nAtoms), nAtoms, [], nConf);
interpMask2 = reshape(repelem(...
    repelem(linspace(1, 0, nConf), 3, 1), 1, nAtoms), nAtoms, [], nConf);
outputCoords = repmat(coords1, 1, 1, nConf) .* interpMask1 + ...
    repmat(coords2, 1, 1, nConf) .* interpMask2;

for i = 1:nConf
    outputStruct.Model(i) = PDBStruct1.Model(1);
    for j = 1:nAtoms
        outputStruct.Model(i).Atom(j).X = outputCoords(j,1,i);
        outputStruct.Model(i).Atom(j).Y = outputCoords(j,2,i);
        outputStruct.Model(i).Atom(j).Z = outputCoords(j,3,i);
    end
end

end

