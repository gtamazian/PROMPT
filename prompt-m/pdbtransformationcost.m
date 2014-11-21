function cost = pdbtransformationcost(PDBStruct, p)
%PDBTRANSFORMATIONCOST Calculates transformation cost.
%   pdbtransformationcost(PDBStruct, p, type) calculates the cost of the 
%   transformation represented by a series of models in the specified 
%   structure PDBStruct. The cost is calculated as a function of the
%   weighted sum of distances between atoms of the adjacent models. The
%   parameter p is the power interatomic distances are raised to.
%
%   The default value for the power p is 2.
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

if nargin < 2
    p = 2;
end

nModels = length(PDBStruct.Model);

% check that the specified structure contains only backbone atoms
for j = 1:nModels
    if ~isempty(setdiff(unique({PDBStruct.Model(j).Atom.AtomName}), ...
            {'N', 'CA', 'C'}))
        error(['the specified PDBStruct object contains models ', ...
            'with non-backbone atoms']);
    end
end

m = atomicmass({PDBStruct.Model(1).Atom.element});  % masses of the atoms

% add masses of side chains and masses of one hydrogen atoms to alpha
% carbons
alphacarbonatoms = PDBStruct.Model(1).Atom( ...
    ismember({PDBStruct.Model(1).Atom.AtomName}, {'CA'}));
m(2:3:end) = m(2:3:end) + sidechainmass({alphacarbonatoms.resName}) + ...
    atomicmass({'H'});

% add masses of hydrogen atoms to nitrogen atoms of the backbone
m(1:3:end) = m(1:3:end) + atomicmass({'H'});

% add masses of oxygen atoms to carbon atoms of the backbone
m(3:3:end) = m(3:3:end) + atomicmass({'O'});

% add a mass of a single hydrogen atom to the N-end of the protein
m(1) = m(1) + atomicmass({'H'});

% add a mass of a hydroxil group to the C-end of the protein
m(end) = m(end) + sum(atomicmass({'O', 'H'}));

% calculate distances between atoms of adjacent models
coords = pdbextractcoords(PDBStruct);
nAtoms = size(coords{1}, 1);
distances = zeros(nAtoms, nModels);    
for j = 2:nModels
    distances(:,j) = sqrt(sum((coords{j} - coords{j-1}).^2, 2));
end
    
cost = sum(m .* sum(distances.^p, 2));

end

