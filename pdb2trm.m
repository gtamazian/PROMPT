function trmodel = pdb2trm(PDBStruct)
%PDB2TRM Create a transformation model from a PDB structure
%   PDB2TRM(PDBStruct) creates a transformation model which represents a
%   transformation from the specified PDB structure.
%
%   See also createmodel trmcreate
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% mail (at) gtamazian (dot) com

coords = pdbextractcoords(PDBStruct);
nModels = size(coords,3);

struct1 = PDBStruct; struct1.Model = struct1.Model(1);
struct2 = PDBStruct; struct2.Model = struct2.Model(end);
trmodel = trmcreate(struct1,struct2,nModels-2);

[~,alpha,psi] = createmodel(PDBStruct);
trmodel.alpha(:,2:end-1) = reduceangles(alpha(:,2:end-1));
trmodel.psi(:,2:end-1) = reduceangles(psi(:,2:end-1));
trmodel = trmupdaterotations(trmodel);

end

