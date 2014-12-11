function h = trmplotfixedrmsd(trmodels,confNum)
%TRMPLOTFIXEDRMSD Plot RMSDs between configurations and the specified one.
%   trmplotfixedrmsd(trmodels,confNum) plots RMSDs between all 
%   configurations of transformations specified in the cell array trmodels 
%   and the configuration with the specified number confNum.
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
rmsdValues = zeros(nTrans,nConf);

for i = 1:nTrans
    coords = trmrestorecoords(trmodels{i});
    for j = 1:nConf
        rmsdValues(i,j) = mean(sqrt(sum((coords{j} - ...
            coords{confNum}).^2,2)));
    end
end

h = plot(transpose(rmsdValues),'-o');
xlabel('Configuration Number');
ylabel('RMSD in AA');

end

