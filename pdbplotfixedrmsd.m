function h = pdbplotfixedrmsd(pdbStruct, modelNo)
%PDBPLOTFIXEDRMSD Plot RMSDs between models and the specified one
%   PDBPLOTFIXEDRMSD(pdbStruct, modelNo) plots RMSDs between all 
%   models of transformations specified in the cell array pdbStruct 
%   and the model with the specified number modelNo.
%
%   See also pdbplotadjrmsd trmplotadjrmsd trmplotfixedrmsd
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014-2016.
% gaik (dot) tamazian (at) gmail (dot) com

% if a single PDB structure is specified instead a cell array, then
% create a cell array with a single element from it
if ~iscell(pdbStruct)
    pdbStruct = {pdbStruct};
end

nTrans = length(pdbStruct);
nModels = length(pdbStruct{1}.Model);
rmsdValues = zeros(nTrans, nModels);

for i = 1:nTrans
    coords = pdbextractcoords(pdbStruct{i});
    for j = 1:nModels
        % superpose the current model to the specified one
        [~, coords{j}] = procrustes(coords{modelNo}, coords{j}, ...
            'scaling', false, 'reflection', false);
        rmsdValues(i,j) = mean(sqrt(sum((coords{j} - ...
            coords{modelNo}).^2,2)));
    end
end

markers = {'s', 'o', '^', '+', 'x', 'd', 'v', '*', 'p', 'h'};
markers = repmat(markers, ceil(nTrans / length(markers)));
markers = markers(:);
for i = 1:size(rmsdValues, 1)
    h = plot(rmsdValues(i,:), ['-', markers{i}]);
    hold on
end
hold off
xlabel('Configuration Number');
ylabel('RMSD in Angstroms');

end

