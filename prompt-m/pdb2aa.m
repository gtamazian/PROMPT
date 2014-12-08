function aa = pdb2aa(PDBStruct)
%PDB2AA Return a cell array of PDB structure amino acid sequences.
%   pdb2aa(PDBStruct) returns a cell array which values are strings of
%   amino acid single letters. Each cell note corresponds to a single model
%   of the specified PDB structure.
%
%   PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

persistent aaabbreviations;
persistent aasingleletter;
persistent aadict;

if any([isempty(aaabbreviations), isempty(aasingleletter), ...
        isempty(aadict)])
    aaabbreviations = {'ALA', 'ARG', 'ASN', 'ASP', 'ASX', 'CYS', 'GLN', ...
        'GLU', 'GLX', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', ...
        'PRO', 'SER', 'THR', 'TRP', 'TYR', 'UNK', 'VAL'};
    aasingleletter = {'A', 'R', 'N', 'D', 'B', 'C', 'Q', 'E', 'Z', 'G', ...
        'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'X', 'V'};
    aadict = containers.Map(aaabbreviations, aasingleletter, ...
        'UniformValues', true);
end

% check if the specified model contains a single model
if length(PDBStruct.Model) > 1
    error('PROMPTM:multModels', 'Multiple models in the PDB structure');
end

% check if specified amino acid names are correct
aminoAcidAbbr = {PDBStruct.Model.Atom.resName};
alphaCarbonIndices = ismember({PDBStruct.Model.Atom.AtomName}, {'CA'});
aminoAcidAbbr = aminoAcidAbbr(alphaCarbonIndices);
missingAminoAcids = setdiff(aminoAcidAbbr, aaabbreviations);
if ~isempty(missingAminoAcids)
    error('PROMPTM:unknownAA', 'Unknown AA abbreviation');
end

% create a cell array of abbreviation cell array
aa = cell2mat(values(aadict, aminoAcidAbbr));

end