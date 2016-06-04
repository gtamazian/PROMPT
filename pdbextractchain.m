function outputPDB = pdbextractchain(PDBStruct, chainID)
%PDBEXTRACTCHAIN Extracts a single chain from a PDB structure
%   PDBEXTRACTCHAIN(PDBStruct, chainID) returns a PDB structure that
%   contains only the specified chain from the provided PDB structure 
%   PDBStruct.
%
%   See also pdbgetchainids
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2016.
% mail (at) gtamazian (dot) com

% check there is a chain with the specified ID
pdbChainIDs = pdbgetchainids(PDBStruct);

if ~ismember(chainID, pdbChainIDs)
    error('PROMPT:pdbextractchain:missingChainID', ...
        'the chain with the specified ID is missing');
end

outputPDB = PDBStruct;
nModels = length(outputPDB.Model);
nAtoms = length(outputPDB.Model(1).Atom);

cellChainID = repmat({chainID}, 1, nAtoms);

for i = 1:nModels
    atomMask = ismember({PDBStruct.Model(1).Atom.chainID}, cellChainID);
    outputPDB.Model(i).Atom = PDBStruct.Model(i).Atom(atomMask);
end

end

