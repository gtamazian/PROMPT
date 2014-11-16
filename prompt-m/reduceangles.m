function Y = reduceangles(X)
%REDUCEANGLES Reduce angular values to [-pi, pi] interval.
%   reduceangles(X) modifies angular values from the given matrix X so that
%   they would lie between -pi and pi.
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

I = X > pi;
while sum(sum(I))
    X(I) = X(I) - 2*pi;
    I = X > pi;
end

I = X < -pi;
while sum(sum(I))
    X(I) = X(I) + 2*pi;
    I = X < -pi;
end

Y = X;

end

