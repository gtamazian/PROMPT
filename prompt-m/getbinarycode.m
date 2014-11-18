function G = getbinarycode(X,N)
%GETBINARYCODE Generates binary numbers from 0 to 2^N.
%   getbinarycode(X,N) returns the matrix which rows are binary numbers
%   from 0 to 2^N. The function should be launched with [] as the X 
%   argument.
%
% PROMPT Toolbox for MATLAB

% By Gaik Tamazian, 2014.
% gaik (dot) tamazian (at) gmail (dot) com

if N > 0
    G = [getbinarycode([X,0], N-1); getbinarycode([X,1], N-1)];
else
    G = X;
end
    
end

