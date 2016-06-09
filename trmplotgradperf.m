function [meanTimes, stdTimes] = trmplotgradperf(model, angleCounts, ...
    nIter, nOutliers, lineSpec)
%TRMPLOTGRADPERF Estimate and plot gradient computation times
%   TRMPLOTGRADPERF(model, angleCounts, nIter, nExtreme, col) produces 
%   a plot showing average gradient computation times for various numbers 
%   of planar and torsion angles specified in angleCounts. The angleCounts 
%   argument is a matrix with two columns: the first one contains the 
%   number of torsion angles and the second one contains the number of 
%   planar angles. nIter specifies the sample size for estimating the 
%   gradient computation time; its default value is 40. The nOutliers
%   argument specifies the number of the least and greatest values from
%   a time sample to be removed as outliers; its default value is 5. The
%   lineSpec argument specifies the plot line style; by default, its value 
%   is '-ob'. The vectors of mean time values and their standard errors
%   are returned.
%
%   See also trmobjfunc trmdistantangleindices
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2016.
% mail (at) gtamazian (dot) com

if nargin < 3
    nIter = 40;
end

if nargin < 4
    nOutliers = 5;
end

if nargin < 5
    lineSpec = '-ob';
end

if 2*nOutliers >= nIter
    error('PROMPT:trmperfplot:tooManyOutliers', ...
        'Too many outliers were specified.');
end

nModels = size(model.r, 2);
nValues = size(angleCounts, 1);
meanTimes = zeros(nValues, 1);
stdTimes = zeros(nValues, 1);

h = waitbar(0, 'Estimating the gradient computation time...');
for i = 1:nValues
    nTorsion = angleCounts(i, 1);
    nPlanar = angleCounts(i, 2);
    torsionIndices = trmdistantangleindices(model, nTorsion, 'torsion');
    planarIndices = trmdistantangleindices(model, nPlanar, 'planar');
    timeValues = zeros(nIter, 1);
    f = @(x) trmobjfunc(model, planarIndices, torsionIndices, x);
    for j = 1:nIter
        angleValues = 2*pi*rand(...
            (nPlanar + nTorsion)*(nModels - 2), 1) - pi;
        tic();
        [~, ~] = f(angleValues);
        timeValues(j) = toc();
        waitbar(((i - 1)*nIter + j)/(nValues * nIter));
    end
    timeValues = sort(timeValues);
    meanTimes(i) = mean(timeValues((nOutliers + 1):(nIter - nOutliers)));
    stdTimes(i) = std(timeValues((nOutliers + 1):(nIter - nOutliers)));
end
close(h);

errors = 1.96 * stdTimes/sqrt(nIter - 2*nOutliers);
errorbar(sum(angleCounts, 2), meanTimes, errors, lineSpec);
xlabel('Number of Angles');
ylabel('Time in Seconds');

end

