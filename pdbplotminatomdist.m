function h = pdbplotminatomdist(pdbStruct)
%PDBPLOTMINATOMDIST Plot minimal interatomic distances.
%   TRMPLOTMINATOMDIST(pdbStruct) plots minimal interatomic distances
%   within models of PDB structures specified in the cell array pdbStruct.
%
%   See also pdbplotadjrmsd pdbplotfixedrmsd trmplotminatomdist
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2016.
% gaik (dot) tamazian (at) gmail (dot) com

% if a single PDB structure model is specified instead a cell array, then
% create a cell array with a single element from it
if ~iscell(pdbStruct)
    pdbStruct = {pdbStruct};
end

nTrans = length(pdbStruct);
nModels = length(pdbStruct{1}.Model);
minDistValues = zeros(nTrans, nModels);

for i = 1:nTrans
    minDistValues(i,:) = pdbmininteratomicdist(pdbStruct{i});
end

markers = {'s', 'o', '^', '+', 'x', 'd', 'v', '*', 'p', 'h'};
markers = repmat(markers, ceil(nTrans / length(markers)));
markers = markers(:);
for i = 1:size(minDistValues, 1)
    h = plot(minDistValues(i,:), ['-', markers{i}]);
    hold on
end
hold off
xlabel('Configuration Number');
ylabel('Minimal Interatomic Distance in Angstroms');

end

