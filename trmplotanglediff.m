function h = trmplotanglediff(trmodel, angleType, sortAngles, ...
    angleNum, markerType)
%TRMPLOTTRANGLEDIFF Plot differences between transformation angles
%   TRMPLOTTRANGLEDIFF(trmodel,angleType,isSorted,angleNum,markerType) 
%   plots absolute values of circular differences between torsion or 
%   planar angles of the first and last configurations of the specified 
%   transformation model trmodel. If sortAngles is true, then the plotted 
%   values are sorted in the descending order. By default, sortAngles is 
%   set to false. The parameter angleNum specifies the number of angles to 
%   be plotted; by default, all angles are shown. The parameter markerType 
%   specifies the plot marker style; its default value is '-'. angleType
%   denotes which angles are considered: planar or torsion. By default,
%   torsion angles are considered.
%
%   See also trmplotadjrmsd trmplotfixedrmsd trmplotminatomdist
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

if nargin < 2
    angleType = 'torsion';
end

if nargin < 3
    sortAngles = false;
end

if nargin < 4
    angleNum = size(trmodel.psi,1);
end

if nargin < 5
    markerType = '-';
end

if strcmp(angleType, 'torsion')
    angleValues = trmodel.psi;
    angleLabel = 'Torsion ';
elseif strcmp(angleType, 'planar')
    angleValues = trmodel.alpha;
    angleLabel = 'Planar ';
else
    error('PROMPT:trmplotanglediff:incorrectAngleType', ...
        'Incorrect angle type specified.');
end

torsAngleDiff = abs(circdist(angleValues(:,1),angleValues(:,end)));
if sortAngles
    torsAngleDiff = sort(torsAngleDiff,'descend');
end
torsAngleDiff = torsAngleDiff(1:angleNum);

h = plot(torsAngleDiff*180/pi,markerType);
axis([0 length(torsAngleDiff) 0 180]);
ylabel('Absolute Circular Distance in Degrees');
if sortAngles
    xlabel([angleLabel, 'Angle Rank']);
else
    xlabel([angleLabel, 'Angle Number']);
end

end

