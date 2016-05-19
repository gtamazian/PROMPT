function [models, optimResults] = trmiteroptim(trmodel, nPlanar, ...
    nTorsion, iterSteps, optimOptions)
%TRMITEROPTIM Iterative optimization scheme.
%   TRMITEROPTIM(trmodel, nPlanar, nTorsion, iterSteps, optimOptions)
%   implements an iterative optimization scheme; the iterSteps array
%   specifies the number of iterations for each stage; nPlanar and nTorsion
%   specify the numbers of planar and torsion angles to be optimized;
%   optimOptions is the structure of the optimization process parameters.
%   The optimization scheme stops if the current stage have not modified
%   the initial point. Also one may specify an empty list as the iterSteps
%   value; in this case, the optimization will continue until local minima
%   are reached at each stage and the current local minimum differs from
%   the ones obtained at the previous stage.
%
%   See also trmobjfunc
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2016.
% mail (at) gtamazian (dot) com

nSteps = length(iterSteps);
models = cell(1, nSteps);
optimResults = cell(1, nSteps);

% first, determine indices of planar and torsion angles to be optimized
P = trmdistantangleindices(trmodel, nPlanar, 'planar');
T = trmdistantangleindices(trmodel, nTorsion, 'torsion');

iStep = 1;
while 1
    currOptimOptions = optimOptions;
    if nSteps > 0
        currOptimOptions.MaxIterations = iterSteps(iStep);
    else
        currOptimOptions.MaxIterations = Inf;
    end
    
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
    
    % check if we should stop the process because the number of steps is
    % equal to the specified one
    if (nSteps > 0) && (iStep == nSteps)
        break
    end
    
    % check if we should stop the process because we have reached the local
    % minimum
    if (iStep > 0) && (output.iterations < 2)
        models = models(1:iStep);
        optimResults = optimResults(1:iStep);
        break
    end
    
    iStep = iStep + 1;
end

