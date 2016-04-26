function minDist = trmmininteratomicdist(trmodel)
%TRMMININTERATOMICDIST Get minimal interatomic distances
%   TRMMININTERATOMICDIST(trmodel) returns the vector of minimal
%   interatomic distances between the atoms within configurations of the
%   specified transformation model. Each vector element corresponds to
%   a certain configuration.
%
%   See also pdbmininteratomicdist
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

coords = trmrestorecoords(trmodel);
minDist = length(coords);

for i = 1:size(coords, 3)
	minDist(i) = min(pdist(coords(:,:,i)));
end

end

