function mass = atomicmass(atomicSymbols)
%ATOMICMASS Get atomic masses
%   ATOMICMASS(atomicSymbols) returns masses of the atoms which symbols 
%   are specified in the cell array atomicSymbols.
%
%   Example:
%
%       % Get atomic masses of carbon, hydrogen and oxygen.
%       atomicmass({'C', 'H', 'O'})
%  
%   See also sidechainmass
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% mail (at) gtamazian (dot) com

persistent symbols;
persistent massValues;
persistent atomicMasses;

if any([isempty(symbols), isempty(massValues), ...
        isempty(atomicMasses)])
    symbols = {'H', 'C', 'N', 'O'};
    massValues = [1.008 12.0107 14.0067 15.9994];
    atomicMasses = containers.Map(symbols, massValues, ...
        'UniformValues', true);
end

mass = cell2mat(values(atomicMasses, atomicSymbols))';

end
