function G = getbinarycode(X, N)
%GETBINARYCODE Generate binary numbers from 0 to 2^N
%   GETBINARYCODE(X, N) returns the matrix which rows are binary numbers
%   from 0 to 2^N. The function should be launched with [] as the X 
%   argument.
%
%   Examples:
%
%       getbinarycode([], 3)
%
%   See also dec2bin
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% mail (at) gtamazian (dot) com

if N > 0
    G = [getbinarycode([X,0], N-1); getbinarycode([X,1], N-1)];
else
    G = X;
end
    
end

