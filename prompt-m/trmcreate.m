function torsmodel = trmcreate(PDBStruct1, PDBStruct2, M)
%TRMCREATE Creates a transformation model.
%   trmcreate(PDBStruct1, PDBStruct2, M) creates a transformation
%   model from two unimodal PDB structures which contain only backbone
%   atoms. The parameter M specifies the number of intermediate
%   configurations in the transformation model.
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

% check if the specified structures contain a single model
if length(PDBStruct1.Model) + length(PDBStruct2.Model) > 2
    error('the specified PDBStruct objects contain more than one model');
end

% remove non-backbone atoms from the specified PDB structures
PDBStruct1 = pdbbackbone(PDBStruct1);
PDBStruct2 = pdbbackbone(PDBStruct2);

% Get atomic masses of backbone atoms.
torsmodel = struct('m', atomicmass({PDBStruct1.Model.Atom.element}), ...
    'StartCoords', atomiccoords(PDBStruct1), ...
    'FinishCoords', atomiccoords(PDBStruct2));

% Process side chains atoms - get a vector of their masses and add them to
% atomic masses of alpha carbons. Also add a mass of one hydrogen atom.
alphacarbonatoms = PDBStruct1.Model.Atom( ...
    ismember({PDBStruct1.Model.Atom.AtomName}, {'CA'}));
torsmodel.m(2:3:end) = torsmodel.m(2:3:end) + ...
    sidechainmass({alphacarbonatoms.resName}) + atomicmass({'H'});

% Add masses of hydrogen atoms to nitrogen atoms of the backbone.
torsmodel.m(1:3:end) = torsmodel.m(1:3:end) + atomicmass({'H'});

% Add masses of oxygen atoms to carbon atoms of the backbone.
torsmodel.m(3:3:end) = torsmodel.m(3:3:end) + atomicmass({'O'});

% Add a mass of a hydrogen atom to N-end of the protein.
torsmodel.m(1) = torsmodel.m(1) + atomicmass({'H'});

% Add a mass of a hydroxil group to C-end of the protein.
torsmodel.m(end) = torsmodel.m(end) + sum(atomicmass({'O', 'H'}));

torsmodel.r     = zeros(length(torsmodel.m) - 1, M+2);
torsmodel.alpha = zeros(length(torsmodel.m) - 2, M+2);
torsmodel.psi   = zeros(length(torsmodel.m) - 3, M+2);

% Calculate bond lengths, planar and torsion angles of the first and last
% configurations.
[torsmodel.r(:,1), torsmodel.alpha(:,1), torsmodel.psi(:,1)] = ...
    createmodel(PDBStruct1);
[torsmodel.r(:,end), torsmodel.alpha(:,end), torsmodel.psi(:,end)] = ...
    createmodel(PDBStruct2);

% Calculate intermediate bond lengths and planar angles.
torsmodel.r(:,2:end-1) = ...
    interpolate(torsmodel.r(:,1), torsmodel.r(:,end), M);
torsmodel.alpha = ...
    circinterp(torsmodel.alpha(:,1), torsmodel.alpha(:,end), M);

% Calculate intermediate torsion angles.
torsmodel.psi = ...
    circinterp(torsmodel.psi(:,1), torsmodel.psi(:,end), M);

% Calculate the rotation matrix for superposition of the second
% conformation to the first one.
torsmodel.U = cell(M+2, 1);
torsmodel.U{1} = eye(3);
prevCoords = torsmodel.StartCoords;
for j = 2:M+2
    currCoords = restorecoords(torsmodel.r(:,j), ...
        torsmodel.alpha(:,j), torsmodel.psi(:,j));
    q = optimquat(prevCoords,currCoords);
    torsmodel.U{j} = quat2rotmat(q);
    r = mean(prevCoords,1) - mean(currCoords*torsmodel.U{j},1);
    prevCoords = currCoords*torsmodel.U{j} + ...
        repmat(r,size(currCoords,1),1);
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

