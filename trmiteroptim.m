function models = trmiteroptim(trmodel, nPlanar, nTorsion, gradTol, ...
    nIterations)
%TRMITEROPTIM Iterative optimization scheme.
%   TRMITEROPTIM(trmodel, nPlanar, nTorsion, schemeOptions, nIterations)
%   implements an iterative optimization scheme; nPlanar and nTorsion
%   specify the numbers of planar and torsion angles to be optimized;
%   gradTol is the gradient threshold value; nIterations is the maximum 
%   number of the optimization process iterations launched at every local
%   step. The function returns a cell array of optimization models; each
%   model in the array corresponds to a stage result.
%
%   See also trmoptim trmobjfunc
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2016.
% mail (at) gtamazian (dot) com

models = cell(1, 1000);

% first, determine indices of planar and torsion angles to be optimized
P = trmdistantangleindices(trmodel, nPlanar, 'planar');
T = trmdistantangleindices(trmodel, nTorsion, 'torsion');

iStep = 0;
iStage = 0;

fprintf('# Stage\t# Step\t# Local step\tFunc value\tGrad norm\n');

while 1
    iStage = iStage + 1;
    jStep = 0;
    exitFlag = 0;
    while ~exitFlag
        iStep = iStep + 1;
        jStep = jStep + 1;
        [x,fval,gnorm,exitFlag] = trmoptim(trmodel, P, T, gradTol, ...
            nIterations);
        trmodel = trmchangeangles(trmodel, P, T, x);
        fprintf('%7d\t%6d\t%12d\t%10e\t%10e\n', iStage, iStep, jStep, ...
            fval, gnorm);
    end
    
    models{iStep} = trmodel;
    
    prevModelCost = trmcost(trmodel);
    trmodel = trmupdaterotations(trmodel);

    if trmcost(trmodel) > prevModelCost
        break;
    end
end

models = models(1:iStep);

end

