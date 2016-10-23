function [x, fval, gnorm, exitFlag] = trmoptim(trmodel, pIndices, ...
    tIndices, gradTol, maxSteps)
%TRMOPTIM Optimize a transformation using the gradient descent
%   [x, fval, gnorm, exitFlag] = trmoptim(trmodel, pIndices, tIndices,
%   gradTol, maxSteps) optimizes the transformation model trmodel by planar
%   and torsion angles which indices are specified by pIndices and
%   tIndices. The optimization stops either when the gradient norm is less
%   than gradTol or the number of optimization steps is greater than
%   maxSteps.
%
%   The function returns the vector of angles x, the corresponding
%   objective function value fval, the corresponding gradient norm value
%   gnorm and the flag exitFlag that indicates why the optimization
%   stopped. If exitFlag is equal to 1, than the gradient norm value was
%   less than gradTol; otherwise, the number of optimization steps exceeded
%   maxSteps.
%
%   See also trmiteroptim trmobjfunc
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2016.
% mail (at) gtamazian (dot) com

    f = @(x) trmobjfunc(trmodel, pIndices, tIndices, x);
    x = trminitialpoint(trmodel, pIndices, tIndices);
    b = 0.5;
    exitFlag = 0;
    
    for iStep = 1:maxSteps
        [fval, gval] = f(x);
        t = 1;
        while (fval - f(x - t*gval)) < (t/2 * norm(gval)^2)
            t = t * b;
        end
        x = x - t*gval;
        if norm(gval) < gradTol
            exitFlag = 1;
            break
        end
    end
    
    [fval, gval] = f(x);
    gnorm = norm(gval);

end

