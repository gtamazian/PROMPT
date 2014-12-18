function h = pdbplotfixedrmsd(pdbStruct, modelNo, fitModels)
%PDBPLOTFIXEDRMSD Plot RMSDs between models and the specified one
%   PDBPLOTFIXEDRMSD(pdbStruct, modelNo) plots RMSDs between all 
%   models of transformations specified in the cell array pdbStruct 
%   and the model with the specified number modelNo. If the
%   parameter fitModels is set true, then RMSDs is calculated for models
%   superposed using the Kabsch transformation.
%
%   See also pdbplotadjrmsd trmplotadjrmsd trmplotfixedrmsd
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

if nargin < 3
    fitModels = false;
end

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
    if fitModels
        for j = 2:nModels
            [~,coords{j}] = procrustes(coords{j-1}, coords{j}, ...
                'scaling', false, 'reflection', false);
        end
    end
    for j = 1:nModels
        rmsdValues(i,j) = mean(sqrt(sum((coords{j} - ...
            coords{modelNo}).^2,2)));
    end
end

h = plot(transpose(rmsdValues),'-o');
xlabel('Configuration Number');
ylabel('RMSD in AA');

end

