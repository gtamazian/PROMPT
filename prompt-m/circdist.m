function d = circdist(X,Y)
%CIRCDIST Return circular distance between two angle matrices.
%   circdist(X,Y), for a pair of matrices with angle values, calculate
%   element-wise circular distances between the angles and return the
%   distance matrix.
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

d = angle(exp(1i*Y)./exp(1i*X));

end
