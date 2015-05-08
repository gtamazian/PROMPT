function aligningAtoms = trmGetAtomsIndsForAlign(trmodel)

nAtoms = size(trmodel.r, 1) + 1;

if ~isfield(trmodel, 'proteinType')
    aligningAtoms = 1:nAtoms;
elseif strcmp(trmodel.proteinType, 'loop')
    aligningAtoms = [1:3, nAtoms-2:nAtoms];
elseif strcmp(trmodel.proteinType, 'n-tail')
    aligningAtoms = nAtoms-2:nAtoms;
elseif strcmp(trmodel.proteinType, 'c-tail')
    aligningAtoms = 1:3;
else
    aligningAtoms = 1:nAtoms;
end

end

