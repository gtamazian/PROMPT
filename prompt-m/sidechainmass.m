function mass = sidechainmass(X)
%AAMASS Return atomic masses of the specified amino acids.
%   aamass(X), for a cell array of amino acid names X, is the vector of
%   their atomic masses.
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

persistent aanotations;
persistent aaweights;
persistent atomicweights;

if any([isempty(aanotations), isempty(aaweights), isempty(atomicweights)])
    aanotations = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', ...
        'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', ...
        'SER', 'THR', 'SEC', 'TRP', 'TYR', 'VAL'};
    aaweights = [71.0779, 156.1857, 114.1026, 115.0874, 103.1429, ...
        129.114, 128.1292, 57.0513, 137.1393, 113.1576, 113.1576, ...
        128.1723, 131.1961, 147.1739, 97.1152, 87.0773, 101.1039, ...
        150.0379, 186.2099, 163.1733, 99.1311];
    atomicweights = containers.Map(aanotations, aaweights, ...
        'UniformValues', true);
end

mass = cell2mat(values(atomicweights, X))';

end
