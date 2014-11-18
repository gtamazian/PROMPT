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

% Calculate bond lengths and planar angles of the initial configuration.
deltax = diff(torsmodel.StartCoords);
torsmodel.r(:,1) = sqrt(sum(deltax.^2, 2));
torsmodel.alpha(:,1) = ...
    acos(dot(deltax(1:end-1,:), deltax(2:end,:), 2) ./ ...
    (torsmodel.r(1:end-1,1) .* torsmodel.r(2:end,1)));
        
% Calculate torsion angles of the initial configuration.
N = cross(deltax(1:end-1,:), deltax(2:end,:), 2);
torsmodel.psi(:,1) = ...
    atan2(torsmodel.r(2:end-1,1) .* ...
    dot(deltax(1:end-2,:), N(2:end,:), 2), ...
    dot(N(1:end-1,:), N(2:end,:), 2));

% Calculate bond lengths and planar angles of the final configuration.
deltax = diff(torsmodel.FinishCoords);
torsmodel.r(:,end) = sqrt(sum(deltax.^2, 2));
torsmodel.alpha(:,end) = ...
    acos(dot(deltax(1:end-1,:), deltax(2:end,:), 2) ./ ...
    (torsmodel.r(1:end-1,end) .* torsmodel.r(2:end,end)));

% Calculate torsion angles of the final configuration.
N = cross(deltax(1:end-1,:), deltax(2:end,:), 2);
torsmodel.psi(:,end) = ...
    atan2(torsmodel.r(2:end-1,end) .* ...
    dot(deltax(1:end-2,:), N(2:end,:), 2), ...
    dot(N(1:end-1,:), N(2:end,:), 2));

% Calculate intermediate bond lengths and planar angles.
torsmodel.r(:,2:end-1) = ...
    interpolate(torsmodel.r(:,1), torsmodel.r(:,end), M);
torsmodel.alpha(:,2:end-1) = ...
    interpolate(torsmodel.alpha(:,1), torsmodel.alpha(:,end), M);

% Calculate intermediate torsion angles.
torsmodel.psi = ...
    circinterp(torsmodel.psi(:,1), torsmodel.psi(:,end), M+2);

% Calculate the rotation matrix for superposition of the second
% conformation to the first one.
torsmodel.U = cell(M+2, 1);
torsmodel.U{1} = eye(3);
prevcoords = restorecoords(torsmodel.r(:,1), torsmodel.alpha(:,1), ...
    torsmodel.psi(:,1));
for j = 2:M+2
    currcoords = restorecoords(torsmodel.r(:,j), ...
        torsmodel.alpha(:,j), torsmodel.psi(:,j));
    [~,~,transformation] = procrustes(prevcoords, currcoords, ....
        'scaling', false, 'reflection', false);
    torsmodel.U{j} = transformation.T;
    prevcoords = currcoords;
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

