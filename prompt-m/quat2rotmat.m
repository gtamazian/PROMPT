function U = quat2rotmat(q)
%QUAT2ROTMAT Convert a quaternion to the rotation matrix.
%   quat2rotmat converts the specified quaternion q to the rotation matrix.
%
%   References:
%       1. Coutsias, Evangelos A., Chaok Seok, and Ken A. Dill. 
%       "Using quaternions to calculate RMSD." Journal of computational 
%       chemistry 25, no. 15 (2004): 1849-1857.
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

U11 = q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2;
U12 = 2*(q(2)*q(3) - q(1)*q(4));
U13 = 2*(q(2)*q(4) + q(1)*q(3));
U21 = 2*(q(2)*q(3) + q(1)*q(4));
U22 = q(1)^2 - q(2)^2 + q(3)^2 - q(4)^2;
U23 = 2*(q(3)*q(4) - q(1)*q(2));
U31 = 2*(q(2)*q(4) - q(1)*q(3));
U32 = 2*(q(3)*q(4) + q(1)*q(2));
U33 = q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2;

U = [U11, U12, U13; U21, U22, U23; U31, U32, U33];

end

