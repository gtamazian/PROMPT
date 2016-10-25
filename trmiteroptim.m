function [models, optimResults] = trmiteroptim(trmodel, nPlanar, ...
    nTorsion, schemeOptions)
%TRMITEROPTIM Iterative optimization scheme.
%   TRMITEROPTIM(trmodel, nPlanar, nTorsion, schemeOptions)
%   implements an iterative optimization scheme; nPlanar and nTorsion
%   specify the numbers of planar and torsion angles to be optimized;
%   schemeOptions is the structure of the optimization process parameters
%   (by default, it is the same as optimoptions('fminunc') with the
%   infinite numbers of iterations and function evaluations);
%   nIterations is the maximum number of optimization stages (the default
%   values is Inf). The optimization scheme stops if the rotation matrix
%   update procedure increases the transformation cost.
%
%   See also trmoptim trmobjfunc
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2016.
% mail (at) gtamazian (dot) com

if nargin < 4
    schemeOptions = optimoptions('fminunc', ...
        'Algorith', 'quasi-newton', ...
        'Display', 'iter', ...
        'OptimalityTolerance', 1e-1, ...
        'SpecifyObjectiveGradient', true, ...
        'MaxIterations', Inf, ...
        'MaxFunctionEvaluations', Inf);
end

models = cell(1, 1000);
optimResults = cell(1, 1000);

% first, determine indices of planar and torsion angles to be optimized
P = trmdistantangleindices(trmodel, nPlanar, 'planar');
T = trmdistantangleindices(trmodel, nTorsion, 'torsion');

iStep = 1;
exitFlag = 0;

while ~exitFlag
    % get the objective function for the current stage
    f = @(x) trmobjfunc(trmodel, P, T, x);

    % get the initial point
    initial_point = trminitialpoint(trmodel, P, T);

    problem.options = schemeOptions;
    problem.solver = 'fminunc';
    problem.x0 = initial_point;
    problem.objective = f;
    
    [x,~,~,output] = fminunc(problem);

    trmodel = trmchangeangles(trmodel, P, T, x);
    models{iStep} = trmodel;
    optimResults{iStep} = output;

    prevModelCost = trmcost(trmodel);
    trmodel = trmupdaterotations(trmodel);
    if (trmcost(trmodel) > prevModelCost) || (output.iterations == 1)
        exitFlag = 1;
    else
        iStep = iStep + 1;
    end
end

models = models(1:iStep);

end

