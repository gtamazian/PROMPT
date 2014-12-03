function trmodel = pdb2trm(PDBStruct)
%PDB2TRM Create a transformation model from a PDB structure object.
%   pdb2trm(PDBStruct) creates a transformation model which represents a
%   transformation from the specified PDB structure.

coords = pdbextractcoords(PDBStruct);
nModels = size(coords,2);

struct1 = PDBStruct; struct1.Model = struct1.Model(1);
struct2 = PDBStruct; struct2.Model = struct2.Model(end);
trmodel = trmcreate(struct1,struct2,nModels-2);

[~,~,psi] = createmodel(PDBStruct);
trmodel.psi(:,2:end-1) = psi(:,2:end-1);

end

