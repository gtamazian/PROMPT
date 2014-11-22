function result = pdbtransformationcompare(PDBStruct1,PDBStruct2)
%PDBTRANSFORMATIONCOMPARE Compare transformations from PDB structures.
%   Detailed explanation goes here

coords1 = pdbextractcoords(PDBStruct1);
coords2 = pdbextractcoords(PDBStruct2);

[r1,alpha1,psi1] = createmodel(PDBStruct1);
[r2,alpha2,psi2] = createmodel(PDBStruct2);

M = length(coords1); % number of configurations in the transformations

% get differences and RMSDs between configuration atom coordinates
coordabsdiff = zeros(M,1);
coordrmsd = zeros(M,1);
rdiff = zeros(M,1);
alphadiff = zeros(M,1);
psidiff = zeros(M,1);
for j = 1:M
    coordabsdiff(j) = norm(coords1{j} - coords2{j});
    coordrmsd(j) = procrustes(coords1{j}, coords2{j}, ...
        'scaling', false, 'reflection', false);
    rdiff(j) = norm(r1 - r2);
    alphadiff(j) = norm(circdist(alpha1,alpha2));
    psidiff(j) = norm(circdist(psi1,psi2));
end

result = struct('coord', coordabsdiff, 'rmsd', coordrmsd, 'r', rdiff, ...
    'alpha', alphadiff, 'psi', psidiff);

end

