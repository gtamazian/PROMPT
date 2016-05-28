function matplot(X, varargin)
%MATPLOT Plot columns of a matrix
%   MATPLOT(X) plots columns of matrix X considering each columns as a
%   separate line in the manner of R's matplot function. A user may
%   customize style of each line being plotted by specifying the following
%   optional arguments.
%
%   - 'type'
%       * 'p' - draw points (default); 
%       * 'l' - draw a line;
%       * 'b' - draw both points and a line;
%   - 'lty' - the line type
%       * '-' - solid (default);
%       * '--' - dashed;
%       * ':' - dotted;
%       * '-.' - dast-dot;
%   - 'lwd' - the line width (0.5 by default);
%   - 'pch' - the marker type
%       * 's' - square (default);
%       * 'o' - circle;
%       * '^' - upward-pointing triangle;
%       * '+' - plus sign;
%       * 'x' - cross;
%       * 'd' - diamond;
%       * 'v' - down-pointing triangle;
%       * '*' - asterisk;
%       * 'p' - pentagram;
%       * 'h' - hexagram;
%   - 'col' - the line and/or marker color
%       * 'k' - black (default);
%       * 'r' - red;
%       * 'g' - green;
%       * 'b' - blue;
%       * 'c' - cyan;
%       * 'm' - magenta;
%       * 'y' - yellow.
%
%   By default, points of varying color and markers are drawn for each
%   column. If a user specifies the style parameter which value is less
%   than the number of columns in matrix X, then these values are coerced
%   to the column number in the same way as in R's matplot function.
%
%   Example
%       matplot(cumsum(rand(10, 5), 2), 'type', 'b', 'col', 'krg', ...
%           'lty', {'-', '--'})

% By Gaik Tamazian, 2016.
% mail (at) gtamazian (dot) com

    typeValues = ['p', 'l', 'b'];
    ltyValues = {'-', '--', ':', '-.'};
    markerValues = ['s', 'o', '^', '+', 'x', 'd', 'v', '*', 'p', 'h'];
    colorValues = ['k', 'r', 'g', 'b', 'c', 'm', 'y'];
    
    p = inputParser;
    p.FunctionName = 'matplot';
    p.PartialMatching = false;
    p.StructExpand = false;
    
    p.addRequired('matrix', @ismatrix);
    p.addParameter('type', 'N', @(x) isempty(setdiff(x, typeValues)));
    p.addParameter('lty', {'N'}, @(x) isempty(setdiff(x, ltyValues)));
    p.addParameter('lwd', 0.5, @iscolumn);
    p.addParameter('pch', 'N', @(x) isempty(setdiff(x, markerValues)));
    p.addParameter('col', 'N', @(x) isempty(setdiff(x, colorValues)));

    parse(p, X, varargin{:});
    
    nLines = size(X, 2);
    
    type = adjustArgument(p.Results.type, 'N', 'p', nLines);
    if ~iscell(p.Results.lty)
        userLty = {p.Results.lty};
    else
        userLty = p.Results.lty;
    end
    lty = adjustArgument(userLty, {'N'}, ltyValues, nLines);
    pch = adjustArgument(p.Results.pch, 'N', markerValues, nLines);
    col = adjustArgument(p.Results.col, 'N', colorValues, nLines);
    lwd = repmat(p.Results.lwd, 1, ceil(nLines / length(p.Results.lwd)));
    
    for iLine = 1:nLines
        switch type(iLine)
            case 'b'
                lineSpec = [lty{iLine}, pch(iLine), col(iLine)];
            case 'l'
                lineSpec = [lty{iLine}, col(iLine)];
            case 'p'
                lineSpec = [pch(iLine), col(iLine)];
        end
        plot(X(:, iLine), lineSpec, 'LineWidth', lwd(iLine));
        hold on
    end
    hold off
    
end

function output = adjustArgument(value, defValue, allValues, argSize)
% This is an auxiliary function that implements the argument coercion for
% the matplot function.

    searchResult = strfind(value, defValue);
    if iscell(searchResult)
        searchResult = cell2mat(searchResult);
    end
    if ~isempty(searchResult)
        value = allValues(1:min(argSize, length(allValues)));         
    end
    output = repmat(value, 1, ceil(argSize / length(value)));
end

