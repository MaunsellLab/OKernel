function [dirName, tableName] = whichData()

% tableName = sprintf('Analysis/Mat Files/masterTable.mat');
% dirName = '/Users/Shared/Data/SCernel/';

tableName = sprintf(' Analysis/Mat Files/masterTable.mat');
if contains(computerName(), 'maunsell')
	dirName = '../../';
else
	dirName = '/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/BehavData/30 PC/';
end

