function backbonePDBStruct = pdbbackbone(PDBStruct, onlyCA)
%PDBBACKBONE Return a PDB structure only with backbone atoms
%   PDBBACKBONE(PDBStruct, onlyCA) returns a PDB structure that contains 
%   only backbone atoms (that is, atoms, which names are N, C, and CA) or
%   only alpha carbon atoms if then onlyCA option is set to true. By
%   default, onlyCA is false.
%
%   See also pdbextractcoords
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014-2016.
% mail (at) gtamazian (dot) com

if nargin < 2
    onlyCA = false;
end

if onlyCA
    atom_names = {'CA'};
else
    atom_names = {'N', 'C', 'CA'};
end

backbonePDBStruct = PDBStruct;
nModels = length(PDBStruct.Model);
for iModel = 1:nModels
    backbonePDBStruct.Model(iModel).Atom = ...
        PDBStruct.Model(iModel).Atom(ismember(...
        {PDBStruct.Model(iModel).Atom.AtomName}, atom_names));
    backbonePDBStruct.Model(iModel).HeterogenAtom = [];
end

end

