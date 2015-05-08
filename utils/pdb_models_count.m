function n = pdb_models_count(pdb_path)
pdb = pdbread(pdb_path);
n = length(pdb.Model);
end