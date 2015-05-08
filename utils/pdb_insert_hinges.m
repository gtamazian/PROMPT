function res = pdb_insert_hinges(orig_pdb, hinge_pdbs, ...
    hinge_torsion_angles_ranges, out_pdb)

pdb = pdbread(orig_pdb);
orig_trm = pdb2trm(pdb,'protein');

for i = 1:length(hinge_torsion_angles_ranges)
    
    h_pdb = pdbread(hinge_pdbs{i});
    h_trm = pdb2trm(h_pdb,'protein');

    orig_trm.psi(hinge_torsion_angles_ranges{i}{1}:...
        hinge_torsion_angles_ranges{i}{end},:) = h_trm.psi;
end

pdb = trm2pdb(orig_trm, pdb);
pdbwrite(out_pdb, pdb);

res = true;

end
