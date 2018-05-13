function trmodel = trmh5load(filename)
%TRMH5LOAD Load a transformation model from the specified HDF5 file.
%   TRMH5LOAD(filename) returns a transformation model loaded from the
%   specified HDF5 file. The following model fields are loaded: 'w' as 'm',
%   'start_coords' as 'StartCoords', 'r' as 'r', 'alpha' as 'alpha',
%   'gamma' as 'psi', and 'rot_mat' as 'U'.
%
%   See also trmh5save pdb2trm
%
% PROMPT Toolbox for MATLAB
% By Gaik Tamazian, 2018.
% mail (at) gtamazian (dot) com

% rotation matrices are written to the HDF5 files transposed by the
% `trmh5save` function - here we transpose them back
tempU = h5read(filename, '/rot_mat');
for k = 1:size(tempU, 3)
    tempU(:,:,k) = transpose(tempU(:,:,k));
end

trmodel = struct('m', h5read(filename, '/w'), ...
    'r', h5read(filename, '/r'), ...
    'alpha', h5read(filename, '/alpha'), ...
    'psi', h5read(filename, '/gamma'), ...
    'StartCoords', transpose(h5read(filename, '/start_coords')), ...
    'U', tempU);

end

