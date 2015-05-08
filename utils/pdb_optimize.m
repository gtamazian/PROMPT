function cost = pdb_optimize(trm_filename, out_pdb_path, ...
    protein_type, n_angles, iter_num, superpose)

pdb = pdbread(trm_filename);
pdb = pdbbackbone(pdb);
pdb_begin = pdb;
pdb_begin.Model = pdb_begin.Model(1);

trm = pdb2trm(pdb, protein_type);

optimized_trm = trm_optimize(trm, n_angles, iter_num, superpose);

cost = trmcost(optimized_trm);
pdb = trm2pdb(optimized_trm, pdb_begin);
pdbwrite(out_pdb_path, pdb);

end

