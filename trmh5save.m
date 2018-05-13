function trmh5save(trmodel, filename)
%TRMH5SAVE Save a transformation model to the specified HDF5 file.
%   TRMH5SAVE(trmodel, filename) writes the transformation model trmodel to
%   the HDF5 file with the name filename. The following model fields are
%   written: 'm' as 'w', 'StartCoords' as 'start_coords', 'r' as 'r',
%   'alpha' as 'alpha', 'psi' as 'gamma', and 'U' and 'rot_mat'.
%
%   See also trmh5load trm2pdb
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2018.
% mail (at) gtamazian (dot) com

if exist(filename, 'file')
    delete(filename);
end

h5create(filename, '/w', size(trmodel.m));
h5write(filename, '/w', trmodel.m);

h5create(filename, '/start_coords', fliplr(size(trmodel.StartCoords)));
h5write(filename, '/start_coords', transpose(trmodel.StartCoords));

h5create(filename, '/r', size(trmodel.r));
h5write(filename, '/r', trmodel.r);

h5create(filename, '/alpha', size(trmodel.alpha));
h5write(filename, '/alpha', trmodel.alpha);

h5create(filename, '/gamma', size(trmodel.psi));
h5write(filename, '/gamma', trmodel.psi);

h5create(filename, '/rot_mat', size(trmodel.U));
h5write(filename, '/rot_mat', trmodel.U);

end

