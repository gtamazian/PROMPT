function U = quat2rotmat(q)
%QUAT2ROTMAT Convert a quaternion to a rotation matrix
%   QUAT2ROTMAT(q) converts the specified quaternion q to a rotation 
%   matrix.
%
%   Example:
%       
%       % Get the quaternion describing the rotation to superpose X and Y
%       % axes to Y and Z axes and convert it to a rotation matrix.
%       rotQuaternion = optimquat([0 0 1; 0 1 0], [0 1 0; 0 0 1])
%       rotMatrix = quat2rotmat(rotQuaternion)
%
%   See also optimquat
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

