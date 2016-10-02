function rmsdValues = pdbplotadjrmsd(pdbStruct, color)
%PDBPLOTADJRMSD Plot RMSDs between adjacent models
%   PDBPLOTADJRMSD(pdbStruct) plots RMSDs between adjacent 
%   models of PDB structures specified in the cell array pdbStruct. A
%   single color for the plot may be specified in the optional color
%   argument.
%
%   See also pdbplotfixedrmsd trmplotadjrmsd trmplotfixedrmsd
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014-2016.
% mail (at) gtamazian (dot) com

% if a single PDB structure model is specified instead a cell array, then
% create a cell array with a single element from it

if nargin < 2
    color = 0;
end

if ~iscell(pdbStruct)
    pdbStruct = {pdbStruct};
end

nTrans = length(pdbStruct);
nModels = length(pdbStruct{1}.Model);
rmsdValues = zeros(nTrans, nModels-1);

for i = 1:nTrans
    coords = pdbextractcoords(pdbStruct{i});
    for j = 2:nModels
        % superpose the current model to the previous one
        [~, coords(:,:,j)] = procrustes(coords(:,:,j-1), coords(:,:,j), ...
            'scaling', false, 'reflection', false);
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
if nModels > 10
    tickStep = 10;
else
    tickStep = nModels - 1;
end
tickPositions = floor(linspace(1, nModels-1, tickStep));
set(ax,'XTick', tickPositions);
xticks = cell(length(tickPositions), 1);
for j = 1:length(tickPositions)
    xticks{j} = [int2str(tickPositions(j)),'-',...
        int2str(tickPositions(j)+1)];
end
set(ax,'XTickLabel',xticks);


end

