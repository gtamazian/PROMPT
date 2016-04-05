function h = trmplotfixedrmsd(trmodels, confNo)
%TRMPLOTFIXEDRMSD Plot RMSDs between configurations and the specified one
%   TRMPLOTFIXEDRMSD(trmodels, confNo) plots RMSDs between all 
%   configurations of transformations specified in the cell array trmodels 
%   and the configuration with the specified number confNo.
%
%   See also trmplotadjrmsd trmplottranglediff trmplotminatomdist
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014-2016.
% gaik (dot) tamazian (at) gmail (dot) com

% if a single transformation model is specified instead a cell array, then
% create a cell array with a single element from it
if ~iscell(trmodels)
    trmodels = {trmodels};
end

nTrans = length(trmodels);
nConf = size(trmodels{1}.r,2);
rmsdValues = zeros(nTrans,nConf);

for i = 1:nTrans
    coords = trmrestorecoords(trmodels{i});
    for j = 1:nConf
        % superpose the current configuration to the specified one
        [~, coords{j}] = procrustes(coords{confNo}, coords{j}, ...
            'scaling', false, 'reflection', false);
        rmsdValues(i,j) = mean(sqrt(sum((coords{j} - ...
            coords{confNo}).^2,2)));
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

