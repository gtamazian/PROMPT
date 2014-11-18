function coordscell = trmrestorecoords(trmodel)
%TRMRESTORECOORDS Restores Cartesian coordinates of transformation atoms.
%   trmrestorecoords(trmodel) returns a cell array of matrices containing
%   Cartesian coordinates of the atoms that constitute the transformation
%   configurations.
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

m = size(trmodel.psi, 2);
n = size(trmodel.r, 1) + 1;
coordscell = cell(1, m);

for i = 1:m
    coordscell{i} = restorecoords(trmodel.r(:,i), ...
        trmodel.alpha(:,i), trmodel.psi(:,i));
    coordscell{i} = coordscell{i}*trmodel.U{i};
    
    % apply the translation
    t = repmat(mean(coordscell{i}, 1), n, 1);
    coordscell{i} = coordscell{i} - t;       
end
