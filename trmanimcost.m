function trmanimcost(trmodel, firstAngle, secondAngle, prefix, ...
    resolution, cycled, caption)
%TRMANIMCOST Procude a series of transformation cost countour plots
%   TRMANIMCOST(model,firstAngle,secondAngle,prefix,resolution,cycled,
%   caption) produces a series of PNG files that show countour plots of the
%   trmodel transformation cost as a function of the specified pair of
%   torsion angles (firstAngle and secondAngle). The output file names are
%   composed of the specified prefix and the number. One may specify the
%   number of points between -pi and pi the countour plot will be build on
%   in the 'resolution' parameter; its default value is 10. Also one may
%   produce a series of figures that correspond to a cycled animation by
%   specifying the 'cycled' parameter equal to true. The plot caption can
%   be  added by specifying its value in the 'caption' parameter.
%
%   See also pdbanimrmch
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2016.
% mail (at) gtamazian (dot) com

if nargin < 5
    resolution = 10;
end

if nargin < 6
    cycled = true;
end

if nargin < 7
    caption = '';
end

nConf = size(trmodel.r, 2);

x = linspace(-pi, pi, resolution);
y = linspace(-pi, pi, resolution);

[X, Y] = meshgrid(x, y);
Z = zeros([size(X), nConf]);

% produce a series of contour plots
point = trminitialpoint(trmodel, [], [firstAngle, secondAngle]);
for k = 2:(nConf - 1)
    for i = 1:resolution
        for j = 1:resolution
            tempPoint = point;
            tempPoint(2*(k - 2) + 1) = X(i, j);
            tempPoint(2*(k - 2) + 2) = Y(i, j);
            Z(i, j, k) = trmobjfunc(trmodel, [], ...
                [firstAngle, secondAngle], tempPoint);
        end
    end
end

% write them to PNG files
for i = 2:(nConf - 1)
    contour(X, Y, Z(:, :, i));
    xlabel(['Torsion Angle #', num2str(firstAngle)]);
    ylabel(['Torsion Angle #', num2str(secondAngle)]);
    title([caption, ' Configuration #', num2str(i)]);
    grid
    hold on
    plot(point(2*(i - 2) + 1), point(2*(i - 2) + 2), 'xr', ...
        'MarkerSize', 10, 'LineWidth', 2);
    hold off
    print('-dpng', [prefix, 'A', num2str(i, '%.2d'), '.png']);
end

if cycled
    for i = (nConf - 1):-1:2
        contour(X, Y, Z(:, :, i));
        xlabel(['Torsion Angle #', num2str(firstAngle)]);
        ylabel(['Torsion Angle #', num2str(secondAngle)]);
        title([caption, ' Configuration #', num2str(i)]);
        grid
        hold on
        plot(point(2*(i - 2) + 1), point(2*(i - 2) + 2), 'xr', ...
            'MarkerSize', 10, 'LineWidth', 2);
        hold off
        print('-dpng', [prefix, 'B', num2str(nConf - i + 1, '%.2d'), ...
            '.png']);
    end
end

end

