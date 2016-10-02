function outputPDB = pdbextractmodel(PDBStruct, modelNo)
%PDBEXTRACTMODEL Extract the specified model from a PDB structure.
%   PDBEXTRACTMODEL(PDBStruct, modelNo) extracts the model with the
%   specified number modelNo from the provided PDB structure PDBStruct.
%
%   See also pdbextractchain pdbextractaltloc
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2016.
% mail (at) gtamazian (dot) com

nModels = length(PDBStruct.Model);
if modelNo > nModels
    error('PROMPT:pdbextractmodel:incorrectModelNumber', ...
        'the model with the specified number is missing');
end

outputPDB = PDBStruct;
outputPDB.Model = PDBStruct.Model(modelNo);

end

