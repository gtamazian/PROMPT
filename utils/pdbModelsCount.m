function n = pdbModelsCount(pdbPath)
pdb = pdbread(pdbPath);
n = length(pdb.Model);
end