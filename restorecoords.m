function coords = restorecoords(r, alpha, psi)
%RESTORECOORDS Restore Cartesian coordinates of atoms
%   restorecoords(r, alpha, psi) restores Cartesian coordinates of the 
%   atoms from the specified bond length, planar and torsion angle vectors
%   (r, alpha and psi vectors, respectively).
%
%   The first point is zero, the second point lies on OX axis and the third
%   one lies on XOY coordinate plane. The Natural Extension Reference frame
%   method is used for placing atoms [1].
%
%   Example:
%
%       r = [1.5044, 1.5330, 1.3314, 1.4580];
%       alpha = [1.2062, 1.1223, 1.0117];
%       psi = [-2.4857, 3.1121];
%       coords = restorecoords(r, alpha, psi)
%
%   References:
%
%       [1] Parsons, J., Holmes, J. B., Rojas, J. M., Tsai, J., & Strauss, 
%       C. E. M. (2005). Practical conversion from torsion space to 
%       Cartesian space for in silico protein synthesis. Journal of %
%       computational chemistry, 26(10), 1063 8. doi:10.1002/jcc.20237
%
%   See also createmodel trmrestorecoords
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% mail (at) gtamazian (dot) com

nAtoms = length(r) + 1;
coords = zeros(nAtoms, 3);

coords(1,:) = [0 0 0];
coords(2,:) = [r(1) 0 0];
coords(3,:) = coords(2,:) + [r(2)*cos(alpha(1)) r(2)*sin(alpha(1)) 0 ];

coords(4:end,1) = r(3:end).*cos(alpha(2:end));
coords(4:end,2) = r(3:end).*sin(alpha(2:end)).*cos(psi);
coords(4:end,3) = r(3:end).*sin(alpha(2:end)).*sin(psi);

for i = 4:nAtoms
    bc = coords(i-1,:) - coords(i-2,:); bc = bc./norm(bc);
    n = cross3d(coords(i-2,:) - coords(i-3,:), bc); n = n./norm(n);
    M = [bc', cross3d(n, bc)', n']; 
    coords(i,:) = (M*coords(i,:)')' + coords(i-1,:);
end

end

