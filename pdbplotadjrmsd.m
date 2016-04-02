function [h, rmsdValues] = pdbplotadjrmsd(pdbStruct)
%PDBPLOTADJRMSD Plot RMSDs between adjacent models
%   PDBPLOTADJRMSD(pdbStruct) plots RMSDs between adjacent 
%   models of PDB structures specified in the cell array pdbStruct.
%
%   See also pdbplotfixedrmsd trmplotadjrmsd trmplotfixedrmsd
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014-2016.
% gaik (dot) tamazian (at) gmail (dot) com

% if a single PDB structure model is specified instead a cell array, then
% create a cell array with a single element from it
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
        [~, coords{j}] = procrustes(coords{j-1}, coords{j}, ...
            'scaling', false, 'reflection', false);
        rmsdValues(i,j-1) = mean(sqrt(sum((coords{j} - ...
            coords{j-1}).^2,2)));
    end
end

markers = {'s', 'o', '^', '+', 'x', 'd', 'v', '*', 'p', 'h'};
for i = 1:size(rmsdValues, 1)
    h = plot(rmsdValues(i,:), ['-', markers{i}]);
    hold on
end
hold off
xlabel('Configuration Pair');
ylabel('RMSD in Angstroms');

% modify x axis tick labels
ax = gca;
xticks = get(ax,'XTickLabel');
for j = 1:nModels-1
    xticks{j} = [int2str(j),'-',int2str(j+1)];
end
set(ax,'XTickLabel',xticks);

end

