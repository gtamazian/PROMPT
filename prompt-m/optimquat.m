function q = optimquat(X, Y)
%OPTIMQUAT Quaternion corresponding to an optimal rotation to fit Y to X
%   OPTIMQUAT(X, Y) returns a 4-element vector representing the quaternion
%   which describes the optimal rotation to fit the points in the matrix
%   Y to the points in the matrix X. The implemented algorithm is 
%   described in [1].
%
%   Example:
%       
%       % Get the quaternion describing the rotation to superpose X and Y
%       % axes to Y and Z axes and convert it to a rotation matrix.
%       rotQuaternion = optimquat([0 0 1; 0 1 0], [0 1 0; 0 0 1])
%       rotMatrix = quat2rotmat(rotQuaternion)
%
%   References:
%
%       [1] Coutsias, Evangelos A., Chaok Seok, and Ken A. Dill. 
%       "Using quaternions to calculate RMSD." Journal of computational 
%       chemistry 25, no. 15 (2004): 1849-1857.
%
%   See also procrustes quat2rotmat
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

nRows = size(X, 1);

X = X - repmat(mean(X, 1), nRows, 1);
Y = Y - repmat(mean(Y, 1), nRows, 1);

R = transpose(X)*Y;

% form the matrix F from Coutsias et al., 2004.
F1 = [R(1,1) + R(2,2) + R(3,3), R(2,3) + R(3,2), R(3,1) - R(1,3), ...
    R(1,2) - R(2,1)];
F2 = [R(2,3) - R(3,2), R(1,1) - R(2,2) - R(3,3), R(1,2) + R(2,1), ...
    R(1,3) + R(3,1)];
F3 = [R(3,1) - R(1,3), R(1,2) + R(2,1), -R(1,1) + R(2,2) - R(3,3), ...
    R(2,3) + R(3,2)];
F4 = [R(1,2) - R(2,1), R(1,3) + R(3,1), R(2,3) + R(3,2), ...
    -R(1,1) - R(2,2) + R(3,3)];

F = [F1; F2; F3; F4];

[eigVectors,eigValues] = eig(F);
[~,eigValueRanks] = sort(diag(eigValues), 'descend');

q = eigVectors(:,eigValueRanks(1));
q = q/sum(q.^2);

end

