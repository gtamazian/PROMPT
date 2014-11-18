function T = trmvariateinterpolations(trmodel,N)
%TRMVARIATEINTERPOLATIONS Variate interpolation modes for torsion angles.
%   trmvariateinterpolations(trmodel,N) returns a cell array of 2^N
%   transformations which differ in the direction the most differing angles
%   in the specified transformation were interpolated.
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

T = cell(N,1);
interpmask = getbinarycode([],N);
m = size(trmodel.r, 2); % the number of configurations

% get indices of N torsion angles which differ at most in the first and
% last configurations of the specified transformation
I = trmdistantangleindices(trmodel,N);
I = I(1:N);

for i = 1:2^N
    newmodel = trmodel;
    newmodel.psi(I,:) = circinterp(newmodel.psi(I,1), ...
        newmodel.psi(I,end), m, transpose(interpmask(i,:)));
    newmodel.psi = reduceangles(newmodel.psi);
    T{i} = newmodel;
end

end

