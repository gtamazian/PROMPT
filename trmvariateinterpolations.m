function transformations = trmvariateinterpolations(trmodel, nAngles)
%TRMVARIATEINTERPOLATIONS Vary interpolation modes for torsion angles
%   TRMVARIATEINTERPOLATIONS(trmodel, nAmgles) returns a cell array of 
%   2^nAngles transformations which differ in the direction the most 
%   differing angles in the specified transformation trmodel were 
%   interpolated.
%
%   See also trmcreate
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% mail (at) gtamazian (dot) com

transformations = cell(nAngles,1);
interpmask = 1 - getbinarycode([],nAngles);
nConf = size(trmodel.r, 2);

% get indices of N torsion angles which differ at most in the first and
% last configurations of the specified transformation
I = trmdistantangleindices(trmodel,nAngles);
I = I(1:nAngles);

for i = 1:2^nAngles
    newModel = trmodel;
    newModel.psi(I,:) = circinterp(newModel.psi(I,1), ...
        newModel.psi(I,end), nConf-2, transpose(interpmask(i,:)));
    newModel.psi = reduceangles(newModel.psi);
    transformations{i} = newModel;
end

end

