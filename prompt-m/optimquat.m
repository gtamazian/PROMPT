function q = optimquat(X,Y)
%OPTIMQUAT Quaternion corresponding to an optimal rotation to fit Y to X.
%   optimquat(X,Y) returns a 4-element vector representing the quaternion
%   which describes the optimal rotation to fit Y to X.
%
%   References:
%       1. Coutsias, Evangelos A., Chaok Seok, and Ken A. Dill. 
%       "Using quaternions to calculate RMSD." Journal of computational 
%       chemistry 25, no. 15 (2004): 1849-1857.
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

n = size(X,1);

X = X - repmat(mean(X,1), n, 1);
Y = Y - repmat(mean(Y,1), n, 1);

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

[V,D] = eig(F);
[~,I] = sort(diag(D), 'descend');

q = V(:,I(1));
q = q/sum(q.^2);

end

