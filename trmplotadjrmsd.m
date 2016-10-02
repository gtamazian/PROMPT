function rmsdValues = trmplotadjrmsd(trmodels, color)
%TRMPLOTADJRMSD Plot RMSDs between adjacent configurations
%   TRMPLOTADJRMSD(trmodels,color) plots RMSDs between adjacent
%   configurations of transformation models specified in a cell array
%   trmodels. A single color for the plot may be specified in the optional
%   color argument.
%
%   See also trmplotfixedrmsd trmplottranglediff trmplotminatomdist
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014-2016.
% mail (at) gtamazian (dot) com

% if a single transformation model is specified instead a cell array, then
% create a cell array with a single element from it

if nargin < 2
    color = 0;
end

if ~iscell(trmodels)
    trmodels = {trmodels};
end

nTrans = length(trmodels);
nConf = size(trmodels{1}.r,2);
rmsdValues = zeros(nTrans,nConf-1);

for i = 1:nTrans
    coords = trmrestorecoords(trmodels{i});
    for j = 2:nConf
        % superpose the current configuration to the previos one
        [~, coords(:,:,j)] = procrustes(coords(:,:,j-1), ...
            coords(:,:,j), 'scaling', false, 'reflection', false);
        rmsdValues(i,j-1) = mean(sqrt(sum((coords(:,:,j) - ...
            coords(:,:,j-1)).^2,2)));
    end
end

if color
    if size(rmsdValues, 1) > 10
        warning('PROMPT:trmplotadjrmsd:repetitiveMarkers', ...
            'repetitive markers of the same color will be shown');
    end
    matplot(rmsdValues', 'type', 'b', 'lty', '-', 'col', color);
else
    matplot(rmsdValues', 'type', 'b', 'lty', '-');
end

xlabel('Configuration Pair');
ylabel('RMSD in Angstroms');

% modify x axis tick labels
ax = gca;
if nConf > 10
    tickStep = 10;
else
    tickStep = nConf - 1;
end
tickPositions = floor(linspace(1, nConf-1, tickStep));
set(ax,'XTick', tickPositions);
xticks = cell(length(tickPositions), 1);
for j = 1:length(tickPositions)
    xticks{j} = [int2str(tickPositions(j)),'-',...
        int2str(tickPositions(j)+1)];
end
set(ax,'XTickLabel',xticks);

end

