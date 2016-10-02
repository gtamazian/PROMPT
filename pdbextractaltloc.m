function outputPDB = pdbextractaltloc(PDBStruct, altLocID)
%PDBEXTRACTALTLOC Extract particular alternative locations.
%   PDBEXTRACTALTLOC(PDBStruct, altLocID) chooses atom positions in the
%   PDB structure PDBStruct that correspond to the alternative location
%   specified by altLocID.
%
%   See also pdbextractmodel pdbextractchain
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2016.
% main (at) gtamazian (dot) com

outputPDB = PDBStruct;
nModels = length(PDBStruct.Model);

for i = 1:nModels
    atomMask = ismember({PDBStruct.Model(i).Atom.altLoc}, {'', altLocID});
    outputPDB.Model(i).Atom = outputPDB.Model(i).Atom(atomMask);
end

end

