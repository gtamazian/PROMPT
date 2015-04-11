function h = trmplottranglediff(trmodel, sortAngles, angleNum, markerType)
%TRMPLOTTRANGLEDIFF Plot differences between transformation torsion angles
%   TRMPLOTTRANGLEDIFF(trmodel,isSorted,angleNum,markerType) plots absolute 
%   values of circular differences between torsion angles of the first and 
%   last configurations of the specified transformation model trmodel. If 
%   sortAngles is true, then the plotted values are sorted in the 
%   descending order. By default, sortAngles is set to false. The parameter
%   angleNum specifies the number of angles to be plotted; by default, all 
%   angles are shown. The parameter markerType specifies the plot marker 
%   style; its default value is '-'.
%
%   See also trmplotadjrmsd trmplotfixedrmsd
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

if nargin < 2
    sortAngles = false;
end

if nargin < 3
    angleNum = size(trmodel.psi,1);
end

if nargin < 4
    markerType = '-';
end

torsAngleDiff = abs(circdist(trmodel.psi(:,1),trmodel.psi(:,end)));
if sortAngles
    torsAngleDiff = sort(torsAngleDiff,'descend');
end
torsAngleDiff = torsAngleDiff(1:angleNum);

h = plot(torsAngleDiff*180/pi,markerType);
axis([0 length(torsAngleDiff) 0 180]);
ylabel('Absolute Circular Distance in Degrees');
if sortAngles
    xlabel('Torsion Angle Rank');
else
    xlabel('Torsion Angle Number');
end

end

