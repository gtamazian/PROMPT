function h = trmplotadjrmsd(trmodels)
%TRMPLOTADJRMSD Plot RMSDs between adjacent configurations
%   TRMPLOTADJRMSD(trmodels) plots RMSDs between adjacent configurations of
%   transformation models specified in a cell array trmodels.
%
%   See also trmplotfixedrmsd trmplottranglediff
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

% if a single transformation model is specified instead a cell array, then
% create a cell array with a single element from it
if ~iscell(trmodels)
    trmodels = {trmodels};
end

nTrans = length(trmodels);
nConf = size(trmodels{1}.r,2);
rmsdValues = zeros(nTrans,nConf-1);

for i = 1:nTrans
    coords = trmrestorecoords(trmodels{i});
    for j = 1:nConf-1
        rmsdValues(i,j) = mean(sqrt(sum((coords{j+1} - coords{j}).^2,2)));
    end
end

h = plot(transpose(rmsdValues),'-o');
xlabel('Configuration Pair');
ylabel('RMSD in AA');

% modify x axis tick labels
ax = gca;
xticks = get(ax,'XTickLabel');
for j = 1:nConf-1
    xticks{j} = [int2str(j),'-',int2str(j+1)];
end
set(ax,'XTickLabel',xticks);

end

