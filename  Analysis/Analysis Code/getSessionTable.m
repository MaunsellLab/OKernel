function [U, dataDirName] = getSessionTable(theTable)
%
% Return a subset of the complete session table.  This function serves as an authoritative selector for the subset
% of ssessions used for the various analyses.  
%

  dataDirName = '/Users/Shared/Data/OKernel/';
	tableDataName = [dataDirName ' Analysis/Processed Files.mat'];
  switch theTable
    case {'All', 'all', 'all steps', 'All Steps', 'All steps', 'all ramps', 'All ramps', 'All Ramps'}
      switch theTable
        case {'All', 'all'}
          limits.rampMS = [0, 500];
        case {'all steps', 'All Steps', 'All steps'}
          limits.rampMS = 0;
        case {'all ramps', 'All ramps', 'All Ramps'}
          limits.rampMS = 500;
      end
      limits.animal = {'All'};
      limits.oneDay = [];
      limits.minSessions = 10;                	% require at least n sessions for each animal
      limits.minTrials = 0;
      limits.criterion = 0;
      limits.minDec = 0.1;
      limits.minDPrime = 0.5;
    case {'Example', 'example'}
      limits.rampMS = 0;
      limits.animal = {'902'};
      limits.oneDay = '2019-10-10';
      limits.minSessions = 10;                	% require at least n sessions for each animal
      limits.minTrials = 0;
      limits.criterion = 0;
      limits.minDec = 0.1;
      limits.minDPrime = 0.5;
    otherwise
      fprintf('getSessionTable: unrecognized table type');
      U = [];
      return;
  end
  [U, ~] = getSubset('normal', dataDirName, tableDataName, limits);
end

