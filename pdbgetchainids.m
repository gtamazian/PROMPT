function chainIDs = pdbgetchainids(PDBStruct)
%PDBGETCHAINIDS Get IDs of chains from a PDB structure
%   PDBGETCHAINIDS(PDBStruct) returns a cell array of chain IDs from the
%   specified PDB structure.
%
%   See also pdbextractchain
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2016.
% mail (at) gtamazian (dot) com

chainIDs = unique({PDBStruct.Model(1).Atom.chainID});

end

