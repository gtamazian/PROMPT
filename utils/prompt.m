function cost = prompt(begin_pdb_path, end_pdb_path, ...
    begin_pdb_ind, end_pdb_ind, out_pdb_path, ...
    n_conf, n_angles, init_point_ind, iter_num, ...
    superpose, protein_type)

% nAngels - quantity of angles which are used as parameters in model
% during optimization. They are the top most changeble angles in model.
%
% protein_type - type of model in pdb. There ara some common types: c-tail,
% n-tail, loop, whole protein.

% Sergey Knyazev
% sergey.n.knyazev@gmail.com

pdb_begin = pdbread(begin_pdb_path);
pdb_begin.Model = pdb_begin.Model(begin_pdb_ind);
pdb_begin = pdbbackbone(pdb_begin);

pdb_end = pdbread(end_pdb_path);
pdb_end.Model = pdb_end.Model(end_pdb_ind);
pdb_end = pdbbackbone(pdb_end);

model = trmcreate(pdb_begin, pdb_end, n_conf, ...
    init_point_ind, protein_type);

optimized_model = trm_optimize(model, n_angles, iter_num, superpose);

cost = trmcost(optimized_model);
pdb = trm2pdb(optimized_model, pdb_begin);
pdbwrite(out_pdb_path, pdb);

end

