function D = circdist(X, Y)
%CIRCDIST Calculate circular distance angles
%   CIRCDIST(X, Y) returns element-wise circular distances between angles 
%   specified in the matrices X and Y.
%
%   Example:
%
%       % Get circular distances between the two angle pairs:
%       % (-pi/4, pi/4), (pi, pi/2)
%       circdist([-pi/4, pi], [pi/4, pi/2])
%
%   See also circinterp
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

D = angle(exp(1i*Y)./exp(1i*X));

end
