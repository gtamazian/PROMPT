function h = trmtorsanglediff(trmodel,isSorted)
%TRMTORSANGLEDIFF Plot the difference between given torsion angles.
%   trmtorsanglediff(trmodel,isSorted) plots absolute values of circular
%   differences between torsion angles of the first and last configurations
%   of the specified transformation model. If isSorted is true, then the
%   plotted values are sorted in the descending order. By default, isSorted
%   is set to false.

if nargin < 2
    isSorted = false;
end

torsAngleDiff = abs(circdist(trmodel.psi(:,1),trmodel.psi(:,end)));
if isSorted
    torsAngleDiff = sort(torsAngleDiff,'descend');
end

h = plot(torsAngleDiff*180/pi);
ylabel('Circular Difference Abs Value in Degrees');
if isSorted
    xlabel('Torsion Angle Rank');
else
    xlabel('Torsion Angle Number');
end

end

