function C = cross3d(A, B)
%CROSS3D Simplified evalulation of a vector cross product for 3D vectors.
%   CROSS3D(A, B) returns the vector cross product of three-dimensional
%   vectors A and B. This implementation is faster than the standard cross
%   routine.
%
%   Examples:
%
%       cross3d([1 0 0], [0 1 0])
%
%   See also cross
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

C = [A(2)*B(3) - A(3)*B(2), ...
    A(3)*B(1) - A(1)*B(3), ...
	A(1)*B(2) - A(2)*B(1)];

end

