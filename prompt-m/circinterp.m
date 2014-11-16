function Z = circinterp(X, Y, M, shortArc)
%CIRCINTERP Circular interpolation of angles by the short or long arc.
%   circinterp(X, Y, M, shortArc) returns the matrix of angles that consist
%   of circular interpolations between the given angle values X and Y. M is
%   the number of interpolated values. shortArc is interpreted as a vector
%   of boolean values, TRUE values correspond to short-arc interpolation,
%   FALSE values to long-arc interpolation.
%
% PROMPT toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

n = size(X, 1);

if nargin < 4
    shortArc = true(n, 1);
end

Mu = repmat(1 - shortArc, 1, M);
I = repmat(0:(M-1), n, 1)/(M-1);
Delta = repmat(circdist(X, Y), 1, M);
X = repmat(X, 1, M);

Z = X + (-2*pi*(sign(Delta).*Mu) + Delta).*I;

end

