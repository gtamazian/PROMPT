function n = init_model_max_ind(first_conf_pdb_path, last_conf_pdb_path,...
    first_conf_ind, last_conf_ind, bidir_ang_treshold)

bidir_ang_treshold = double(bidir_ang_treshold);

first_conf = pdbread(first_conf_pdb_path);
first_conf.Model = first_conf.Model(first_conf_ind);

last_conf = pdbread(last_conf_pdb_path);
last_conf.Model = last_conf.Model(last_conf_ind);

trm = trmcreate(first_conf, last_conf, 0);
n = 2^sum(abs(circdist(trm.psi(:,1),trm.psi(:,end)))>=bidir_ang_treshold);
end
