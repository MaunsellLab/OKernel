function [dirName, tableName, projName] = whichData()

% projName = 'JDC';
projName = 'JJC';

tableName = sprintf(' Analysis/Mat Files/masterTable%s.mat', projName);
dirName = '/Users/Shared/Data/OKernel/';
