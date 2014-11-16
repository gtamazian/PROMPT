function [r, alpha, psi] = createmodel(PDBStruct)
%CREATEMODEL Calculate bond lengths, planar and torsion angles.
%   [r, alpha, psi] = createmodel(PDBStruct) calculates bond lengths,
%   planar and torsion angles of the transformation specified by the
%   structure PDBStruct. The bond lengths are returned in the matrix r, the
%   planar angles - in the matrix alpha, and the torsion angles - in the
%   matrix psi. Each column of the matrices r, alpha, and psi corresponds 
%   to the separate configuration of the transformation.
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

nModels = length(PDBStruct.Model);
nAtoms  = length(PDBStruct.Model(1).Atom);

% initialize matrices of bond lengths, planar and torsion angles
r       = zeros(nAtoms-1, nModels);
alpha   = zeros(nAtoms-2, nModels);
psi     = zeros(nAtoms-3, nModels);

% calculate bond lengths and planar and torsion angles for each
% configuration
for i = 1:nModels
    % In matrices r, alpha and psi, the ith column corresponds to the ith
    % configuration.
    deltax = diff([[PDBStruct.Model(i).Atom.X]' ...
        [PDBStruct.Model(i).Atom.Y]' [PDBStruct.Model(i).Atom.Z]']);
    r(:,i) = sqrt(sum(deltax.^2, 2));
    alpha(:,i) = acos(dot(deltax(1:end-1,:), deltax(2:end,:), 2) ./ ...
        (r(1:end-1,i) .* r(2:end,i)));
    N = cross(deltax(1:end-1,:), deltax(2:end,:), 2);
    psi(:,i) = atan2(r(2:end-1,i) .* ...
        dot(deltax(1:end-2,:), N(2:end,:), 2), ...
        dot(N(1:end-1,:), N(2:end,:), 2));
end

end
