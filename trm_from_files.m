function trm = trm_from_files(begin_pdb_path, ...
    end_pdb_path, begin_pdb_ind, end_pdb_ind, n_conf)

begin_pdb = pdbread(begin_pdb_path);
end_pdb = pdbread(end_pdb_path);

first_conf = begin_pdb;
first_conf.Model = first_conf.Model(begin_pdb_ind);
last_conf = end_pdb;
last_conf.Model = last_conf.Model(end_pdb_ind);

trm = trmcreate(first_conf, last_conf, n_conf);

end

