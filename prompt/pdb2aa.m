function aa = pdb2aa(PDBStruct)
%PDB2AA Return a cell array of PDB structure amino acid sequences
%   PDB2AA(PDBStruct) returns a cell array which values are strings of
%   amino acid single letters. Each cell element corresponds to a single
%   model of the specified PDB structure.
%
%   See also sidechainmass
%
%   PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

persistent symbols;
persistent letters;
persistent aminoAcidDict;

if any([isempty(symbols), isempty(letters), isempty(aminoAcidDict)])
    symbols = {'ALA', 'ARG', 'ASN', 'ASP', 'ASX', 'CYS', 'GLN', ...
        'GLU', 'GLX', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', ...
        'PRO', 'SER', 'THR', 'TRP', 'TYR', 'UNK', 'VAL'};
    letters = {'A', 'R', 'N', 'D', 'B', 'C', 'Q', 'E', 'Z', 'G', ...
        'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'X', 'V'};
    aminoAcidDict = containers.Map(symbols, letters, ...
        'UniformValues', true);
end

% check if the specified model contains a single model
if length(PDBStruct.Model) > 1
    error('Multiple models in the PDB structure');
end

% check if specified amino acid names are correct
aminoAcidSymbols = {PDBStruct.Model.Atom.resName};
alphaCarbonIndices = ismember({PDBStruct.Model.Atom.AtomName}, {'CA'});
aminoAcidSymbols = aminoAcidSymbols(alphaCarbonIndices);
missingAminoAcids = setdiff(aminoAcidSymbols, symbols);
if ~isempty(missingAminoAcids)
    error('Unknown AA symbol');
end

% create a cell array of abbreviation cell array
aa = cell2mat(values(aminoAcidDict, aminoAcidSymbols));

end

