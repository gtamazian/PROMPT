function [models, optimResults] = trmiteroptim(trmodel, nPlanar, ...
    nTorsion, iterSteps, optimOptions)
%TRMITEROPTIM Iterative optimization scheme.
%   TRMITEROPTIM(trmodel, nPlanar, nTorsion, iterSteps, optimOptions)
%   implements an iterative optimization scheme; the iterSteps array
%   specifies the number of iterations for each stage; nPlanar and nTorsion
%   specify the numbers of planar and torsion angles to be optimized;
%   optimOptions is the structure of the optimization process parameters.
%
%   See also trmobjfunc
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2016.
% gaik (dot) tamazian (at) gmail (dot) com

nSteps = length(iterSteps);
models = cell(nSteps, 1);
optimResults = cell(nSteps, 1);

% first, determine indices of planar and torsion angles to be optimized
P = trmdistantangleindices(trmodel, nPlanar, 'planar');
T = trmdistantangleindices(trmodel, nTorsion, 'torsion');

for iStep = 1:nSteps
    currOptimOptions = optimOptions;
    currOptimOptions.MaxIterations = iterSteps(iStep);
    
    % get the objective function for the current stage
    f = @(x) trmobjfunc(trmodel, P, T, x);
    
    % get the initial point
    initial_point = trminitialpoint(trmodel, P, T);
    
    [x,~,~,output] = fmincon(f, initial_point, [], [], [], [], ...
        -pi*ones(size(initial_point)), pi*ones(size(initial_point)), ...
        [], currOptimOptions);
    
    trmodel = trmchangeangles(trmodel, P, T, x);
    models{iStep} = trmodel;
    optimResults{iStep} = output;
    
    trmodel = trmupdaterotations(trmodel);
end

