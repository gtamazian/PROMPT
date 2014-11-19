function I = trmdistantangleindices(trmodel, N)
%TRMDISTANTANGLEINDICES Get indices of N most distant torsion angles.
%   trmdistantangleindices(trmodel, N) returns indices of N torsion
%   angles phi and psi which values in the original conformations differ at
%   most.
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

[~,I] = sort(abs(circdist(trmodel.psi(:,1), trmodel.psi(:,end))), ...
    'descend');
I = I(1:N);

end

