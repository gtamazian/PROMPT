function mindist = trmmininteratomicdist(trmodel)
%TRMMININTERATOMICDIST Returns the minimal interatomic distances.
%   trmmininteratomicdist(trmodel) returns the vector of the minimal
%   interatomic distances between the atoms that belong to the same
%   transformation configuration. Each element of the vector corresponds to
%   the certain configuration.
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

coords = trmrestorecoords(trmodel);
mindist = length(coords);

for i = 1:length(coords)
	mindist(i) = min(pdist(coords{i}));
end

end

