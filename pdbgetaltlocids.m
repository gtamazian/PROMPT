function altLocIDs = pdbgetaltlocids(PDBStruct)
%PDBGETALTLOCIDS Get alternative loci IDs from a PDB structure
%   PDBGETALTLOCIDS(PDBStruct) returns a cell array of alternative loci IDs
%   from the specified PDB structure PDBStruct.
%
%   See also pdbextractaltloc
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2016.
% mail (at) gtamazian (dot) com

altLocIDs = unique({PDBStruct.Model(1).Atom.altLoc});

end

