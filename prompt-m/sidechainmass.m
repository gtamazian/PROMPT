function result = sidechainmass(sideChains)
%SIDECHAINMASS Get atomic masses of protein side chains
%   SIDECHAINMASS(sideChains) returns a vector of total atomic masses for
%   side chains which names are specified by the cell array sideChains.
%
%   Example:
%
%       sidechainmass({'MET', 'GLY', 'CYS'})
%
%   See also atomicmass
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

persistent symbols;
persistent masses;
persistent sideChainMasses;

if any([isempty(symbols), isempty(masses), isempty(sideChainMasses)])
    symbols = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', ...
        'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', ...
        'SER', 'THR', 'SEC', 'TRP', 'TYR', 'VAL'};
    masses = [71.0779, 156.1857, 114.1026, 115.0874, 103.1429, ...
        129.114, 128.1292, 57.0513, 137.1393, 113.1576, 113.1576, ...
        128.1723, 131.1961, 147.1739, 97.1152, 87.0773, 101.1039, ...
        150.0379, 186.2099, 163.1733, 99.1311];
    sideChainMasses = containers.Map(symbols, masses, ...
        'UniformValues', true);
end

result = cell2mat(values(sideChainMasses, sideChains))';

end
