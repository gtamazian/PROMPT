function [F, G] = trmobjfunc(trmodel, angleIndices, angles)
%TRMOBJFUNC Calculate the cost function and its gradient for the model
%   TRMOBJFUNC(trmodel,indices,angles) returns the cost of the specified
%   transformation trmodel for the values of torsion angles specified in
%   the vector angles. The angles are also defined by their indices in the
%   given vector angleIndices.
%
%   See also trmcost
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

trmodel.psi(angleIndices,2:end-1) = ...
    reshape(angles, length(angleIndices), size(trmodel.psi, 2) - 2);

[F, X] = trmcost(trmodel);

if nargout > 1
    G = g(trmodel, X);
    G = G(angleIndices, 2:end-1);
end

end

function G = g(model,X)
%G Calculate the transformation cost gradient vector.
%   g(model,X) returns the gradient vector of the transformation. It uses a
%   precalculated matrix of Cartesian coordinates X of transformation atoms
%   for better performance.
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

D = d(X);

n = size(model.r, 1) + 1;
m = size(model.r, 2);

G = zeros(n-3, m);

for j = 2:m-1
    H = 2*X{j} - X{j+1} - X{j-1};
    for i = 1:n-3
        G(i,j) = sum(2*model.m .* sum(H .* D{j,i}, 2));
    end
end

end

function D = d(X)
%D Calculate values for a transformation cost gradient vector.
%   d(X) calculates values for a transformation cost gradient vector based
%   on a matrix of transformation atom Cartesian coordinates.
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

m = length(X); % the number of configurations
n = size(X{1}, 1); % the number of atoms

U = u(X);
D = cell(m,n-3);

mask = [zeros(3,n-3); tril(ones(n-3,n-3))];

for j = 1:m
    for i = 1:n-3
        D{j,i} = cross3dvec(repmat(U{j}(i,:), n, 1), ...
            X{j} - repmat(X{j}(i+1,:), n, 1));
        D{j,i} = D{j,i} .* repmat(mask(:,i), 1, 3);
        D{j,i} = D{j,i} - repmat(sum( ...
            cross3dvec(repmat(U{j}(i,:), n-i-2, 1), X{j}(i+3:end,:) - ...
            repmat(X{j}(i+1,:), n-i-2, 1)), 1)/n, n, 1);
    end
end

end

function U = u(X)
%U Calculate values for a transformation cost gradient vector.
%   u(X) calculates values for a transformation cost gradient vector based
%   on a matrix of transformation atom Cartesian coordinates.
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

m = length(X); % the number of configurations
n = size(X{1}, 1); % the number of atoms

U = cell(1, m);
for i = 1:m
    U{i} = X{i}(3:n-1,:) - X{i}(2:n-2,:);
    normU = repmat(sqrt(sum(U{i}.^2, 2)), 1, 3);
    U{i} = U{i} ./ normU;
end

end

function c = cross3dvec(a,b)
%CROSS3DVEC Calculate cross products for matrix rows.
%   cross3dvec(a,b) returns a matrix of cross products for rows of the
%   specified matrices.
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

c = [a(:,2).*b(:,3)-a(:,3).*b(:,2)
    a(:,3).*b(:,1)-a(:,1).*b(:,3)
    a(:,1).*b(:,2)-a(:,2).*b(:,1)];

c = reshape(c, size(a));

end

