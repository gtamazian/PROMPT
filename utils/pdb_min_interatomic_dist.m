function d = pdb_min_interatomic_dist(pdb_file_name, protein_type)
pdb = pdbread(pdb_file_name);
trm = pdb2trm(pdb, protein_type);
d = min(trmmininteratomicdist(trm));
end

