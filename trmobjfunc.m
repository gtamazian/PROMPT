function [F, G] = trmobjfunc(trmodel, planarIndices, torsionIndices, ...
    angles)
%TRMOBJFUNC Calculate the cost function and its gradient for the model
%   TRMOBJFUNC(trmodel,planarIndices,torsionIndices,angles) returns the 
%   cost of the specified transformation trmodel for the values of torsion
%   angles specified in the vector angles. The angles are also defined by 
%   their indices in the vector arguments planarIndices and torsionIndices.
%
%   See also trmcost
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014-2016.
% gaik (dot) tamazian (at) gmail (dot) com

trmodel = trmchangeangles(trmodel, planarIndices, torsionIndices, angles);

[F, xTrans, xConf] = trmcost(trmodel);

if nargout > 1
    G = g(trmodel, planarIndices, torsionIndices, xTrans, xConf);
end

end

function G = g(model, planarIndices, torsionIndices, xTrans, xConf)
%G Calculate the transformation cost gradient vector.
%   g(model, planarIndices, torsionIndices, xTrans, xConf) returns the 
%   gradient vector of the transformation.
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014-2016.
% gaik (dot) tamazian (at) gmail (dot) com

nAtoms = size(xTrans, 1);
nConf = size(xTrans, 3);
G = zeros(length(planarIndices) + length(torsionIndices), nConf-2);

vS = s(xTrans);
vR = r(xConf);
vN = n(vR);
vP = p(vR);

maskP = reshape(repelem(tril(ones(nAtoms, nAtoms-1), -2), 1, 3), ...
    [nAtoms, 3, nAtoms-1]);
maskT = reshape(repelem(tril(ones(nAtoms, nAtoms-1), -3), 1, 3), ...
    [nAtoms, 3, nAtoms-1]);
    
for j = 2:nConf-1
    vQ = q(xConf(:,:,j));
        
    vQP = vQ .* maskP;
    vQT = vQ .* maskT;
    
    for angleNumber = 1:length(planarIndices)
        i = planarIndices(angleNumber);
        G(angleNumber,j-1) = 2 * sum(model.m .* dot(vS(:,:,j), ...
            cross(repmat(vP(i,:,j), nAtoms, 1), vQP(:,:,i) - ...
            repmat(mean(vQP(:,:,i)), [nAtoms, 1])) * model.U(:,:,j), 2));
    end
    
    for angleNumber = 1:length(torsionIndices)
        i = torsionIndices(angleNumber);
        G(length(planarIndices) + angleNumber,j-1) =  ...
            2 * sum(model.m .* dot(vS(:,:,j), ...
            cross(repmat(vN(i,:,j), nAtoms, 1), vQT(:,:,i) - ...
            repmat(mean(vQT(:,:,i)), [nAtoms, 1])) * model.U(:,:,j), 2));
    end
end

end

function output = s(x)
 
    output = zeros(size(x));
    output(:,:,2:end-1) = 2*x(:,:,2:end-1) - x(:,:,1:end-2) - x(:,:,3:end);

end

function output = r(x)

    output = x(2:end,:,:) - x(1:end-1,:,:);

end

function output = n(r)

    output = r ./ repmat(sqrt(sum(r.^2, 2)), [1, 3, 1]);
    output = output(2:end,:,:);

end

function output = p(r)

    tempP = cross(r(1:end-1,:,:), r(2:end,:,:), 2);
    tempPNorm = repmat(sqrt(sum(tempP.^2, 2)), [1, 3, 1]);
    output = tempP ./ tempPNorm;

end

function output = q(x)

    nAtoms = size(x, 1);
    output = repmat(x, [1, 1, nAtoms-1]) - ...
        permute(reshape(repelem(x(2:end,:), nAtoms, 1), ...
        [nAtoms, nAtoms-1, 3]), [1, 3, 2]);
    
end


