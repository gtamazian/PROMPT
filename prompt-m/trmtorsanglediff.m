function h = trmtorsanglediff(trmodel,isSorted,angleNum,markerType)
%TRMTORSANGLEDIFF Plot the difference between given torsion angles.
%   trmtorsanglediff(trmodel,isSorted,angleNum,markerType) plots absolute 
%   values of circular differences between torsion angles of the first and 
%   last configurations of the specified transformation model. If isSorted 
%   is true, then the plotted values are sorted in the descending order. 
%   By default, isSorted is set to false. The parameter angleNum specifies 
%   the number of angles to be plotted; by default, all angles are shown.
%   The parameter markerType specifies the plot marker style; its default
%   value is '-'.
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

if nargin < 2
    isSorted = false;
end

if nargin < 3
    angleNum = size(trmodel.psi,1);
end

if nargin < 4
    markerType = '-';
end

torsAngleDiff = abs(circdist(trmodel.psi(:,1),trmodel.psi(:,end)));
if isSorted
    torsAngleDiff = sort(torsAngleDiff,'descend');
end
torsAngleDiff = torsAngleDiff(1:angleNum);

h = plot(torsAngleDiff*180/pi,markerType);
ylabel('Circular Difference Abs Value in Degrees');
if isSorted
    xlabel('Torsion Angle Rank');
else
    xlabel('Torsion Angle Number');
end

end

