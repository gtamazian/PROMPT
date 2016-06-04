function minDistValues = pdbplotminatomdist(pdbStruct, color)
%PDBPLOTMINATOMDIST Plot minimal interatomic distances.
%   TRMPLOTMINATOMDIST(pdbStruct) plots minimal interatomic distances
%   within models of PDB structures specified in the cell array pdbStruct.
%   A single color for the plot may be specified in the optional color
%   argument.
%
%   See also pdbplotadjrmsd pdbplotfixedrmsd trmplotminatomdist
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2016.
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
minDistValues = zeros(nTrans, nModels);

for i = 1:nTrans
    minDistValues(i,:) = pdbmininteratomicdist(pdbStruct{i});
end

if color
    if size(minDistValues, 1) > 10
        warning('PROMPT:trmplotadjrmsd:repetitiveMarkers', ...
            'repetitive markers of the same color will be shown');
    end
    matplot(minDistValues', 'type', 'b', 'lty', '-', 'col', color);
else
    matplot(minDistValues', 'type', 'b', 'lty', '-');
end
xlabel('Configuration Number');
ylabel('Minimal Interatomic Distance in Angstroms');

end

