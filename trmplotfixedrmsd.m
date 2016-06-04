function rmsdValues = trmplotfixedrmsd(trmodels, confNo, color)
%TRMPLOTFIXEDRMSD Plot RMSDs between configurations and the specified one
%   TRMPLOTFIXEDRMSD(trmodels, confNo, color) plots RMSDs between all
%   configurations of transformations specified in the cell array trmodels 
%   and the configuration with the specified number confNo. A single color
%   for the plot may be specified in the optional color argument.
%
%   See also trmplotadjrmsd trmplottranglediff trmplotminatomdist
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014-2016.
% mail (at) gtamazian (dot) com

% if a single transformation model is specified instead a cell array, then
% create a cell array with a single element from it

if nargin < 3
    color = 0;
end

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
        [~, coords(:,:,j)] = procrustes(coords(:,:,confNo), coords(:,:,j), ...
            'scaling', false, 'reflection', false);
        rmsdValues(i,j) = mean(sqrt(sum((coords(:,:,j) - ...
            coords(:,:,confNo)).^2,2)));
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

