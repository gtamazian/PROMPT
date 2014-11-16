function PDBStruct = trm2pdb(trmodel, initialPDBStruct)
%TRM2PDB Converts a transformation model to a PDB struct. 
%   trm2pdb(trmodel, initialPDBStruct) returns the PDB structure that
%   correposonds to the transformation presented in the transformation
%   model trmodel. The initial PDB structure the transformartion has been 
%   created from is also required.
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

PDBStruct = initialPDBStruct;
PDBStruct.Model = PDBStruct.Model(1);
nModels = size(trmodel.psi, 2);
coordscell = trmrestorecoords(trmodel);

for i = 1:nModels
    PDBStruct.Model(i) = initialPDBStruct.Model(1);
    for j = 1:size(coordscell{i}(:,1), 1)
        PDBStruct.Model(i).Atom(j).X = coordscell{i}(j,1);
        PDBStruct.Model(i).Atom(j).Y = coordscell{i}(j,2);
        PDBStruct.Model(i).Atom(j).Z = coordscell{i}(j,3);
    end
    PDBStruct.Model(i).MDLSerNo = i;
end

end

