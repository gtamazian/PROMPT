function reducedAngles = reduceangles(angles)
%REDUCEANGLES Reduce angular values to the interval [-pi, pi].
%   REDUCEANGLES(angles) modifies angular values from the specified matrix
%   angles so that they would lie between -pi and pi.
%
%   Example:
%
%       reduceangles([3*pi/2, 2*pi, 0])
%
%   See also circdist
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

I = angles > pi; % indices of angles to be reduced
while sum(sum(I))
    angles(I) = angles(I) - 2*pi;
    I = angles > pi;
end

I = angles < -pi;
while sum(sum(I))
    angles(I) = angles(I) + 2*pi;
    I = angles < -pi;
end

reducedAngles = angles;

end

