function minDistValues = trmplotminatomdist(trmodels, color)
%TRMPLOTMINATOMDIST Plot minimal interatomic distances.
%   TRMPLOTMINATOMDIST(trmodels) plots minimal interatomic distances
%   within configurations of the specified transformations. A single color
%   for the plot may be specified in the optional color argument.
%
%   See also trmplotadjrmsd trmplotfixedrmsd pdbplotminatomdist
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2016.
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
minDistValues = zeros(nTrans,nConf);

for i = 1:nTrans
    minDistValues(i,:) = trmmininteratomicdist(trmodels{i});
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
