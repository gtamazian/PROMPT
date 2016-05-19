function Z = circinterp(X, Y, nValues, shortArc)
%CIRCINTERP Circular interpolation of angles by the specified arcs
%   CIRCINTERP(X, Y, nValues, shortArc) returns the matrix that contains
%   circular interpolations between the specified angle vectors X and Y.
%   nValue is the number of interpolated values. shortArc is a boolean
%   vector which specifies directions of circular interpolation for
%   each pair of angles from X and Y. True values correspond to short-arc
%   interpolation, false values to long-arc interpolation.
%
%   Examples:
%
%       % Interpolate three angles by the short arc.
%       circinterp([pi/4, pi/2, pi], [-pi/4, pi, 0], 4)
%
%       % Interpolate the first pair of angles by the short arc and the
%       % second one by the long arc.
%       circinterp([pi/4, pi/2], [-pi/4, pi], 4, [true, false])
%
%   See also circdist
%
% PROMPT toolbox for MATLAB

% By Gaik Tamazian, 2014.
% mail (at) gtamazian (dot) com

nAngles = numel(X);

if nargin < 4
    shortArc = true(nAngles, 1);
end

% if X, Y or shortArc are row vectors, convert them to column vectors
if size(X, 1) == 1
    X = transpose(X);
end
if size(Y, 1) == 1
    Y = transpose(Y);
end
if size(shortArc, 1) == 1
    shortArc = transpose(shortArc);
end

% check if X and Y are column vectors and have the same number of elements
xIsVector = size(X, 2) == 1;
yIsVector = size(Y, 2) == 1;
sameNumberElements = size(X, 1) == size(Y, 1);
if ~(xIsVector && yIsVector && sameNumberElements)
    error('Incorrect angle vectors specified.')
end

Mu = repmat(1 - shortArc, 1, nValues+2);
I = repmat(0:(nValues+1), nAngles, 1)/(nValues+1);
circDistances = repmat(circdist(X, Y), 1, nValues+2);
X = repmat(X, 1, nValues+2);

Z = X + (-2*pi*(sign(circDistances).*Mu) + circDistances).*I;

end

