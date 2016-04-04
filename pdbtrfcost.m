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

% check atoms in the provided PDB structures
atomNames = {PDBStruct.Model(1).Atom.AtomName};

if all(strcmp(unique(atomNames), {'CA'}))
    onlyCA = true;
elseif all(strcmp(unique(atomNames), {'C', 'CA', 'N'}))
    onlyCA = false;
else
    error('PROMPT:trmcreate:nonBackboneAtoms', ...
        'Non-backbone atoms in the provided PDB structures.');
end

atomicMasses = atomicmass({PDBStruct.Model(1).Atom.element});

% add masses of side chains and masses of one hydrogen atoms to alpha
% carbons
alphaCarbonAtoms = PDBStruct.Model(1).Atom( ...
    ismember({PDBStruct.Model(1).Atom.AtomName}, {'CA'}));

if onlyCA
    atomicMasses = atomicMasses + ...
        sidechainmass({alphaCarbonAtoms.resName}) + ...
        2*atomicmass({'H'}) + sum(atomicmass({'N', 'C', 'O'}));
else
    % Process side chains atoms - get a vector of their masses and add
    % them to atomic masses of alpha carbons. Also add a mass of one 
    % hydrogen atom.
    atomicMasses(2:3:end) = atomicMasses(2:3:end) + ...
        sidechainmass({alphaCarbonAtoms.resName}) + atomicmass({'H'});

    % Add masses of hydrogen atoms to nitrogen atoms of the backbone.
    atomicMasses(1:3:end) = atomicMasses(1:3:end) + atomicmass({'H'});

    % Add masses of oxygen atoms to carbon atoms of the backbone.
    atomicMasses(3:3:end) = atomicMasses(3:3:end) + atomicmass({'O'});
end

% Add a mass of two hydrogen atoms to N-end of the protein.
atomicMasses(1) = atomicMasses(1) + 2*atomicmass({'H'});

% Add a mass of a hydroxil group to C-end of the protein.
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

