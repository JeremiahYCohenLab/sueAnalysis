function [root, sep] = currComputer(customRoot)

if nargin > 0 && ~isempty(customRoot)
    root = customRoot;
else 
    root = 'F:\';
end

% Default behavior based on OS
if ismac
    sep = '/';
elseif ispc
    sep = '\';
else
    error('Unsupported OS');
end
