function h = trmplotminatomdist(trmodels)
%TRMPLOTMINATOMDIST Plot minimal interatomic distances.
%   TRMPLOTMINATOMDIST(trmodels) plots minimal interatomic distances
%   within configurations of the specified transformations.
%
%   See also trmplotadjrmsd trmplotfixedrmsd pdbplotminatomdist
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2016.
% mail (at) gtamazian (dot) com

% if a single transformation model is specified instead a cell array, then
% create a cell array with a single element from it
if ~iscell(trmodels)
    trmodels = {trmodels};
end

nTrans = length(trmodels);
nConf = size(trmodels{1}.r,2);
minDistValues = zeros(nTrans,nConf);

for i = 1:nTrans
    minDistValues(i,:) = trmmininteratomicdist(trmodels{i});
end

matplot(minDistValues', 'type', 'b', 'lty', '-');
xlabel('Configuration Number');
ylabel('Minimal Interatomic Distance in Angstroms');

end
