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

if length(PDBStruct.Model) > 1
    error('the specified PDBStruct object contains multiple models')
end

backbonePDBStruct = PDBStruct;
backbonePDBStruct.Model.Atom = PDBStruct.Model.Atom( ...
    ismember({PDBStruct.Model.Atom.AtomName}, {'N' 'C' 'CA'}));
backbonePDBStruct.Model.HeterogenAtom = [];

end

