function [d,Z,A] = sievefit(X,Y,nMinPoints,epsDistance,epsRotation)
%SIEVEFIT Apply the sieve-fit procedure to fit Y to X.
%   sievefit(X,Y,nMinPoints,eps) performs the sieve-fit procedure to fit 
%   points from the matrix Y to the points from the matrix X. Three 
%   stopping criteria are provided: 
%       1. by the number of atoms (nMinPoints),
%       2. by the difference between subsets being fitted (epsDistance),
%       3. by the norm of difference between rotation matrices at the 
%       current and previous steps (epsRotation).
%   If any stopping criterion is satisfied, the sieve-fit procedure is
%   stopped. 
%
%   The function returns the distance d between the fitted points (d), the 
%   matrix Z of transformed points from Y, and the matrix A describing the 
%   obtained rotation.
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

nAtoms = size(X,1);

% check the specified stopping criteria
if nargin < 4
    nMinPoints = round(nAtoms/2);
end

if nargin < 5
    epsDistance = 0;
end

if nargin < 6
    epsRotation = 0;
end

tempX = X;
tempY = Y;

prevRotation = eye(1);
while nAtoms >= nMinPoints
    [d,Z,transformation] = procrustes(tempX,tempY,'scaling',false, ...
        'reflection',false);
    
    if d < epsDistance
        break
    end
    
    if norm(transformation.T - prevRotation) < epsRotation
        break
    end
    
    % determine the pair of points which are the most distant from each
    % other and remove them
    distances = sum((tempX - Z).^2,2);
    [~,I] = max(distances);
    indices = setdiff(1:nAtoms,I);
    tempX = tempX(indices,:);
    tempY = tempY(indices,:);
    nAtoms = nAtoms - 1;
    
    prevRotation = transformation.T;
end

A = transformation.T;
r = mean(X*transformation.T,1) - mean(Y,1);
Z = X*A - repmat(r,size(X,1),1);
d = sum(sqrt(sum((X-Z).^2,2)));

end

