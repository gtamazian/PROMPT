function h = pdbplotadjrmsd(pdbStruct, fitModels)
%PDBPLOTADJRMSD Plot RMSDs between adjacent models
%   PDBPLOTADJRMSD(pdbStruct, fitModels) plots RMSDs between adjacent 
%   models of PDB structures specified in the cell array pdbStruct. If the
%   parameter fitModels is set true, then RMSDs is calculated for models
%   superposed using the Kabsch transformation.
%
%   See also pdbplotfixedrmsd trmplotadjrmsd trmplotfixedrmsd
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

if nargin < 2
    fitModels = false;
end

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
    if fitModels
        for j = 2:nModels
            [~,coords{j}] = procrustes(coords{j-1}, coords{j}, ...
                'scaling', false, 'reflection', false);
        end
    end
    for j = 1:nModels-1
        rmsdValues(i,j) = mean(sqrt(sum((coords{j+1} - coords{j}).^2,2)));
    end
end

h = plot(transpose(rmsdValues),'-o');
xlabel('Configuration Pair');
ylabel('RMSD in AA');

% modify x axis tick labels
ax = gca;
xticks = get(ax,'XTickLabel');
for j = 1:nModels-1
    xticks{j} = [int2str(j),'-',int2str(j+1)];
end
set(ax,'XTickLabel',xticks);

end

