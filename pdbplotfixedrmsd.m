function rmsdValues = pdbplotfixedrmsd(pdbStruct, modelNo, color)
%PDBPLOTFIXEDRMSD Plot RMSDs between models and the specified one
%   PDBPLOTFIXEDRMSD(pdbStruct, modelNo) plots RMSDs between all 
%   models of transformations specified in the cell array pdbStruct 
%   and the model with the specified number modelNo. A single color for
%   the plot may be specified in the optional color argument.
%
%   See also pdbplotadjrmsd trmplotadjrmsd trmplotfixedrmsd
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014-2016.
% mail (at) gtamazian (dot) com

% if a single PDB structure is specified instead a cell array, then
% create a cell array with a single element from it

if nargin < 3
    color = 0;
end

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
        [~, coords(:,:,j)] = procrustes(coords(:,:,modelNo), ...
            coords(:,:,j), 'scaling', false, 'reflection', false);
        rmsdValues(i,j) = mean(sqrt(sum((coords(:,:,j) - ...
            coords(:,:,modelNo)).^2,2)));
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
xlabel('Configuration Number');
ylabel('RMSD in Angstroms');

end

