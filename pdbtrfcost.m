function cost = pdbtrfcost(PDBStruct)
%PDBTRFCOST Calculate transformation cost from a PDB structure
%   PDBTRFCOST(PDBStruct) calculates the cost of the transformation 
%   represented by a series of models in the specified PDB structure 
%   PDBStruct. The cost is calculated as a function of the weighted sum of 
%   squared distances between atoms of the adjacent models.
%
%   See also trmcost
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

p = 2;
nModels = length(PDBStruct.Model);

% check that the specified structure contains only backbone atoms
for iModel = 1:nModels
    if ~isempty(...
            setdiff(unique({PDBStruct.Model(iModel).Atom.AtomName}), ...
            {'N', 'CA', 'C'}))
        error(['the specified PDBStruct object contains models ', ...
            'with non-backbone atoms']);
    end
end

atomicMasses = atomicmass({PDBStruct.Model(1).Atom.element});

% add masses of side chains and masses of one hydrogen atoms to alpha
% carbons
alphacarbonatoms = PDBStruct.Model(1).Atom( ...
    ismember({PDBStruct.Model(1).Atom.AtomName}, {'CA'}));
atomicMasses(2:3:end) = atomicMasses(2:3:end) + ...
    sidechainmass({alphacarbonatoms.resName}) + ...
    atomicmass({'H'});

% add masses of hydrogen atoms to nitrogen atoms of the backbone
atomicMasses(1:3:end) = atomicMasses(1:3:end) + atomicmass({'H'});

% add masses of oxygen atoms to carbon atoms of the backbone
atomicMasses(3:3:end) = atomicMasses(3:3:end) + atomicmass({'O'});

% add a mass of a single hydrogen atom to the N-end of the protein
atomicMasses(1) = atomicMasses(1) + atomicmass({'H'});

% add a mass of a hydroxil group to the C-end of the protein
atomicMasses(end) = atomicMasses(end) + sum(atomicmass({'O', 'H'}));

% calculate distances between atoms of adjacent models
coords = pdbextractcoords(PDBStruct);
nAtoms = size(coords{1}, 1);
atomicDist = zeros(nAtoms, nModels);    
for iModel = 2:nModels
    atomicDist(:,iModel) = sqrt(sum((coords{iModel} - ...
        coords{iModel-1}).^2, 2));
end
    
cost = sum(atomicMasses .* sum(atomicDist.^p, 2));

end

