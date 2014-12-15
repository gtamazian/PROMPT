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
% gaik (dot) tamazian (at) gmail (dot) com

persistent symbols;
persistent massValues;
persistent atomicMasses;

if any([isempty(symbols), isempty(massValues), ...
        isempty(atomicMasses)])
    symbols = {'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', ...
        'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', ...
        'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', ...
        'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', ...
        'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', ...
        'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La'};
    massValues = [1 4 6.94 9.01 10.81 12.01 14.01 16 19 20.18 23 24.31 ...
        26.98 28.09 30.97 32.07 35.45 39.95 39.1 40.08 44.96 47.88 ...
        50.94 52 54.94 55.85 58.93 58.69 63.55 65.39 69.72 72.59 74.92 ...
        78.96 79.9 83.8 85.47 87.62 88.91 91.22 92.91 95.94 97.91 ...
        101.07 102.91 106.42 107.87 112.41 114.82 118.71 121.75 127.6 ...
        126.91 131.29 132.91 137.33 138.96];
    atomicMasses = containers.Map(symbols, massValues, ...
        'UniformValues', true);
end

mass = cell2mat(values(atomicMasses, atomicSymbols))';

end
