function cost = pdbtransformationcost(PDBStruct, p)
%PDBTRANSFORMATIONCOST Calculates transformation cost.
%   pdbtransformationcost(PDBStruct, p, type) calculates the cost of the 
%   transformation represented by a series of models in the specified 
%   structure PDBStruct. The cost is calculated as a function of the
%   weighted sum of distances between atoms of the adjacent models. The
%   parameter p is the power interatomic distances are raised to.
%
%   There are two cost types. In the type A cost, the sum of interatomic
%   distances is raised to the power p. In the type B cost, interatomic
%   distances are raised to the power p and then summed.
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
for i = 1:nModels
    if ~isempty(setdiff(unique({PDBStruct.Model(i).Atom.AtomName}), ...
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
coordsprev = [[PDBStruct.Model(1).Atom.X]' ...
        [PDBStruct.Model(1).Atom.Y]' [PDBStruct.Model(1).Atom.Z]'];
nAtoms = size(coordsprev, 1);
distances = zeros(nAtoms, nModels);    
for i = 2:nModels
    coordscur = [[PDBStruct.Model(i).Atom.X]' ...
        [PDBStruct.Model(i).Atom.Y]' [PDBStruct.Model(i).Atom.Z]'];
    distances(:,i) = sqrt(sum((coordscur - coordsprev).^2, 2));
    coordsprev = coordscur;
end
    
cost = sum(m .* sum(distances, 2).^p);

end

