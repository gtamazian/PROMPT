function enPDBStruct = pdbenrich(origPDBStruct, fullPDBStruct)
%PDBENRICH Enrich a PDB structure with fields from another one.
%   PDBENRICH(origPDBStruct,fullPDBStruct) adds to the origPDBStruct
%   structure fields rhat are missing in it but present in the
%   fullPDBStruct structure.
%
%   See also pdbstrip
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2016.
% mail (at) gtamazian (dot) com

commonFields = intersect(fieldnames(fullPDBStruct), ...
    fieldnames(origPDBStruct));
for i = 1:length(commonFields)
    fullPDBStruct.(commonFields{i}) = origPDBStruct.(commonFields{i});
end

enPDBStruct = fullPDBStruct;

end

