function d = pdbMinInteratomicDist(pdbFileName)
pdb = pdbread(pdbFileName);
trm = pdb2trm(pdb);
d = min(trmmininteratomicdist(trm));
end

