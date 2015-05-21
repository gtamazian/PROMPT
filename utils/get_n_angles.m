function n_angles = get_n_angles(begin_pdb_path, ...
    end_pdb_path, begin_pdb_ind, end_pdb_ind, angles_threshold)

trmodel = trm_from_files(begin_pdb_path, end_pdb_path, ...
    begin_pdb_ind, end_pdb_ind, 0);

n_angles = length(find(abs(circdist(trmodel.psi(:,1), ...
    trmodel.psi(:,end))) > angles_threshold));

end

