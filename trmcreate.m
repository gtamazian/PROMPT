function trmodel = trmcreate(PDBStruct1, PDBStruct2, nConf)
%TRMCREATE Create a transformation model
%   TRMCREATE(PDBStruct1, PDBStruct2, nModels) creates a transformation
%   model from two PDB structures which contain a single model and
%   correspond to the same protein. Note that amino acid content of both 
%   structures must be the same. The parameter nConf specifies the number 
%   of intermediate configurations in the transoformation model.
%
%   See also createmodel pdb2trm
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

% check if the specified structures contain a single model
if length(PDBStruct1.Model) + length(PDBStruct2.Model) > 2
    error('the specified PDBStruct objects contain more than one model');
end

% check atoms in the provided PDB structures
atomNames1 = {PDBStruct1.Model.Atom.AtomName};
atomNames2 = {PDBStruct2.Model.Atom.AtomName};

if length(atomNames1) ~= length(atomNames2)
    error('PROMPT:trmcreate:differentAtomNumbers', ...
        'Different atom numbers in the provided PDB structures.');
end

if ~all(strcmp(atomNames1, atomNames2))
    error('PROMPT:trmcreate:differentAtomNames', ...
        'Different atom names in the provided PDB structures.');
end

if all(strcmp(unique(atomNames1), {'CA'}))
    onlyCA = true;
elseif all(strcmp(unique(atomNames1), {'C', 'CA', 'N'}))
    onlyCA = false;
else
    error('PROMPT:trmcreate:nonBackboneAtoms', ...
        'Non-backbone atoms in the provided PDB structures.');
end

% Get atomic masses of the atoms.
trmodel = struct('m', atomicmass({PDBStruct1.Model.Atom.element}), ...
    'StartCoords', atomiccoords(PDBStruct1), ...
    'FinishCoords', atomiccoords(PDBStruct2));

alphaCarbonAtoms = PDBStruct1.Model.Atom( ...
    ismember({PDBStruct1.Model.Atom.AtomName}, {'CA'}));

if onlyCA
    trmodel.m = trmodel.m + sidechainmass({alphaCarbonAtoms.resName}) + ...
        2*atomicmass({'H'}) + sum(atomicmass({'N', 'C', 'O'}));
else
    % Process side chains atoms - get a vector of their masses and add
    % them to atomic masses of alpha carbons. Also add a mass of one 
    % hydrogen atom.
    trmodel.m(2:3:end) = trmodel.m(2:3:end) + ...
        sidechainmass({alphaCarbonAtoms.resName}) + atomicmass({'H'});

    % Add masses of hydrogen atoms to nitrogen atoms of the backbone.
    trmodel.m(1:3:end) = trmodel.m(1:3:end) + atomicmass({'H'});

    % Add masses of oxygen atoms to carbon atoms of the backbone.
    trmodel.m(3:3:end) = trmodel.m(3:3:end) + atomicmass({'O'});
end

% Add a mass of two hydrogen atoms to N-end of the protein.
trmodel.m(1) = trmodel.m(1) + 2*atomicmass({'H'});

% Add a mass of a hydroxil group to C-end of the protein.
trmodel.m(end) = trmodel.m(end) + sum(atomicmass({'O', 'H'}));

trmodel.r     = zeros(length(trmodel.m) - 1, nConf+2);
trmodel.alpha = zeros(length(trmodel.m) - 2, nConf+2);
trmodel.psi   = zeros(length(trmodel.m) - 3, nConf+2);

% Calculate bond lengths, planar and torsion angles of the first and last
% configurations.
[trmodel.r(:,1), trmodel.alpha(:,1), trmodel.psi(:,1)] = ...
    createmodel(PDBStruct1);
[trmodel.r(:,end), trmodel.alpha(:,end), trmodel.psi(:,end)] = ...
    createmodel(PDBStruct2);

% Calculate intermediate bond lengths and planar angles.
trmodel.r(:,2:end-1) = ...
    interpolate(trmodel.r(:,1), trmodel.r(:,end), nConf);
trmodel.alpha = ...
    circinterp(trmodel.alpha(:,1), trmodel.alpha(:,end), nConf);

% Calculate intermediate torsion angles.
trmodel.psi = ...
    circinterp(trmodel.psi(:,1), trmodel.psi(:,end), nConf);

% Calculate the rotation matrix for superposition of the second
% conformation to the first one.
trmodel.U = zeros(3, 3, nConf + 2);
trmodel.U(:,:,1) = eye(3);
prevCoords = trmodel.StartCoords;
for j = 2:nConf+2
    currCoords = restorecoords(trmodel.r(:,j), ...
        trmodel.alpha(:,j), trmodel.psi(:,j));
    [~,prevCoords,t] = procrustes(prevCoords, currCoords, ...
        'scaling', false, 'reflection', false);
    trmodel.U(:,:,j) = t.T;
end

end

function result = interpolate(start, finish, M)
    result = ...
        repmat((M:-1:1)/(M+1), length(start), 1).*repmat(start, 1, M) + ...
        repmat((1:M)/(M+1), length(finish), 1).*repmat(finish, 1, M);
end

function coords = atomiccoords(PDBStruct)
    coords = [[PDBStruct.Model.Atom.X]' [PDBStruct.Model.Atom.Y]' ...
        [PDBStruct.Model.Atom.Z]'];
end

