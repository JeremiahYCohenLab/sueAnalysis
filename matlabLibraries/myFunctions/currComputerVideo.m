function [root, sep] = currComputerVideo()

if ismac
    root = '/Volumes/cooper/';
  %  root = '/Volumes/bbari1/';
    sep = '/';
elseif ispc
   root = 'D:\video\';
  %  root = 'Z:\';
%    root = 'C:\Users\zhixi\Documents\data\';
  %  root = 'D:\';
%     root = 'C:\Users\zhixi\Documents\data\';
    sep = '\';
end