function backbonePDBStruct = pdbbackbone(PDBStruct)
%PDBBACKBONE Return a PDB structure only with backbone atoms
%   PDBBACKBONE(PDBStruct) returns a PDB structure that contains only
%   backbone atoms (that is, atoms, which names are N, C, and CA).
%
%   See also pdbextractcoords
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

backbonePDBStruct = PDBStruct;
nModels = length(PDBStruct.Model);
for iModel = 1:nModels
    backbonePDBStruct.Model(iModel).Atom = ...
        PDBStruct.Model(iModel).Atom(ismember(...
        {PDBStruct.Model(iModel).Atom.AtomName}, {'N' 'C' 'CA'}));
    backbonePDBStruct.Model(iModel).HeterogenAtom = [];
end

end

