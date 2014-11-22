function mindist = pdbmininteratomicdist(PDBStruct)
%PDBMININTERATOMICDIST Returns the minimal interatomic distances.
%   pdbmininteratomicdist(PDBStruct) returns the vector of the minimal
%	interatomic distances between the atoms that belong to the same
%	transformation configuration. Each element of the vector corresponds
%	to the certain configuration.
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

coords = pdbextractcoords(PDBStruct);
mindist = length(coords);

for i = 1:length(coords)
	mindist(i) = min(pdist(coords{i}));
end

end

