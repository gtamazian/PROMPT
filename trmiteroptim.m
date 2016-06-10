function [models, optimResults] = trmiteroptim(trmodel, nPlanar, ...
    nTorsion, optimOptions, nIterations)
%TRMITEROPTIM Iterative optimization scheme.
%   TRMITEROPTIM(trmodel, nPlanar, nTorsion, optimOptions, nIterations)
%   implements an iterative optimization scheme; nPlanar and nTorsion
%   specify the numbers of planar and torsion angles to be optimized;
%   optimOptions is the structure of the optimization process parameters;
%   nIterations is the maximum number of optimization stages (the default
%   values is Inf). The optimization scheme stops if the current stage
%   have not modified the initial point or the number of stages is greater
%   tha  the nIterations value.
%
%   See also trmobjfunc
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2016.
% mail (at) gtamazian (dot) com

if nargin < 5
    nIterations = Inf;
end

models = cell(1, 1000);
optimResults = cell(1, 1000);

% first, determine indices of planar and torsion angles to be optimized
P = trmdistantangleindices(trmodel, nPlanar, 'planar');
T = trmdistantangleindices(trmodel, nTorsion, 'torsion');

iStep = 1;
while 1
    % get the objective function for the current stage
    f = @(x) trmobjfunc(trmodel, P, T, x);
    
    % get the initial point
    initial_point = trminitialpoint(trmodel, P, T);
    
    problem.options = optimOptions;
    problem.solver = 'fminunc';
    problem.x0 = initial_point;
    problem.objective = f;

    [x,~,~,output] = fminunc(problem);
    
    trmodel = trmchangeangles(trmodel, P, T, x);
    models{iStep} = trmodel;
    optimResults{iStep} = output;
    
    trmodel = trmupdaterotations(trmodel);

    % check if we should stop the process because we have reached the local
    % minimum
    if ((iStep > 0) && (output.iterations < 2)) || (iStep == nIterations)
        models = models(1:iStep);
        optimResults = optimResults(1:iStep);
        break
    end
    
    iStep = iStep + 1;
end

