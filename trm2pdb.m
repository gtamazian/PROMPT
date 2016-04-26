function PDBStruct = trm2pdb(trmodel, initialPDBStruct)
%TRM2PDB Convert a transformation model to a PDB structure object
%   TRM2PDB(trmodel, initialPDBStruct) returns the PDB structure that
%   corresponds to the transformation presented in the specified 
%   transformation model trmodel. The initial PDB structure the 
%   transformartion has been created from is also required.
%
%   See also trmcreate
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

PDBStruct = initialPDBStruct;
PDBStruct.Model = PDBStruct.Model(1);
nModels = size(trmodel.psi, 2);
coords = trmrestorecoords(trmodel);

for i = 1:nModels
    PDBStruct.Model(i) = initialPDBStruct.Model(1);
    for j = 1:size(coords(:,1,i), 1)
        PDBStruct.Model(i).Atom(j).X = coords(j,1,i);
        PDBStruct.Model(i).Atom(j).Y = coords(j,2,i);
        PDBStruct.Model(i).Atom(j).Z = coords(j,3,i);
    end
end

% add the model numbers; we do it in a separate loop to be able to assign
% model structures from initialPDBStruct to PDBStruct in the loop above
for i = 1:nModels
    PDBStruct.Model(i).MDLSerNo = i;
end

end

