function strPDBStruct = pdbstrip(PDBStruct)
%PDBSTRIP Strip a PDB structure of all fields except for atom coordinates.
%   PDBSTRIP(PDBStruct) returns a PDB structure derived from the specified
%   one by removing all its fields except atom coordinates, namely,
%   Model.Atoms. If the atom coordinate field is absent, then the error
%   message is shown.
%
%   See also pdbenrich
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2016.
% mail (at) gtamazian (dot) com

if isfield(PDBStruct, 'Model') && isfield(PDBStruct.Model, 'Atom')
    strPDBStruct = struct;
    strPDBStruct.Model.Atom = PDBStruct.Model.Atom;
else
    error('PROMPT:pdbstrip:missingAtomCoordinates', ...
        'Atom coordinates missing in the provided PDB structure.');

end

