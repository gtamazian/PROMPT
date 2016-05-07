function [coords, confCoords] = trmrestorecoords(trmodel)
%TRMRESTORECOORDS Restore Cartesian coordinates of transformation atoms
%   TRMRESTORECOORDS(trmodel) returns two arrays of matrices of atom
%   Cartesian coordinates: coords contains coorditanes of the
%   transformation atoms and confCoords contains coordinates of the 
%   unaligned configuration atoms. Both arrays store configurations in
%   the last (third) dimension.
%
%   See also restorecoords
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014-2016.
% gaik (dot) tamazian (at) gmail (dot) com

nConf = size(trmodel.psi, 2);
nAtoms = size(trmodel.r, 1) + 1;
coords = zeros(nAtoms, 3, nConf);
confCoords = zeros(nAtoms, 3, nConf);
coords(:,:,1) = trmodel.StartCoords;
firstConfTranslation = repmat(mean(coords(:,:,1), 1), nAtoms, 1);

for i = 2:nConf
    confCoords(:,:,i) = restorecoords(trmodel.r(:,i), ...
        trmodel.alpha(:,i), trmodel.psi(:,i));
    coords(:,:,i) = confCoords(:,:,i)*trmodel.U(:,:,i);
    
    % apply the translation
    currTranslation = repmat(mean(coords(:,:,i), 1), nAtoms, 1);
    coords(:,:,i) = coords(:,:,i) - currTranslation + firstConfTranslation;       
end

end

