function pdbanimrmch(PDBStruct, prefix, cycled)
%PDBANIMRMCH Produce a series of Ramachandran plots for PDB models
%   PDBANIMRMCH(PDBStruct,prefix,cycled) produces a series of PNG files
%   that show the Ramachandran plots for models of the specified PDB
%   structure 'PDBStructure'. The output file names are composed of the
%   specified 'prefix' and the number. Using the 'cycled' parameter, one
%   may produce a series of figures corresponding a cycled animation by
%   specifying it equal to true.
%
%   See also trmanimcost
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2016.
% mail (at) gtamazian (dot) com

if nargin < 3
    cycled = true;
end

nModels = length(PDBStruct.Model);

% the header is required for the ramachandran function
if ~isfield(PDBStruct, 'Header') || ~isfield(PDBStruct.Header, 'idCode')
    PDBStruct.Header.idCode = '';
end

for iModel = 1:nModels
    tempStruct = PDBStruct;
    tempStruct.Model = PDBStruct.Model(iModel);
    ramachandran(tempStruct);
    title(['Configuration #', num2str(iModel)]);
    print('-dpng', [prefix, 'A', num2str(iModel, '%.2d'), '.png']);
end

if cycled
    for iModel = nModels:-1:1
        tempStruct = PDBStruct;
        tempStruct.Model = PDBStruct.Model(iModel);
        ramachandran(tempStruct);
        title(['Configuration #', num2str(iModel)]);
        print('-dpng', [prefix, 'B', num2str(nModels - iModel + 1, ...
            '%.2d'), '.png']);
    end
end

end

