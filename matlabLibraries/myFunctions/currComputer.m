function [root, sep] = currComputer()

if ismac
    root = '/Volumes/cooper/';
  %  root = '/Volumes/bbari1/';
    sep = '/';
elseif ispc
<<<<<<< HEAD
   root = 'D:\';
  %  root = 'Z:\';
%    root = 'C:\Users\zhixi\Documents\data\';
=======
  %  root = 'D:\';
    root = 'Y:\';
%     root = 'C:\Users\zhixi\Documents\data\';
>>>>>>> cfe5d4b5f3b3d76ac28bea90f69606a89622482a
    sep = '\';
end