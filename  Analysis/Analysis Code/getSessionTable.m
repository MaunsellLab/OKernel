function [U, dataDirName, limits] = getSessionTable(theSubset)
%
% Return a subset of the complete session table.  This function serves as an authoritative selector for the subset
% of ssessions used for the various analyses.  
%
% unfiltered -- no selection criteria (except numStim > 0)
% all -- standard selection criteria

  dataDirName = '/Users/Shared/Data/OKernel/';
  load([dataDirName ' Analysis/Mat Files/masterTable.mat'], 'T');
  limits = setLimits(theSubset);
  if isempty(limits)
    T = [];
    return;
  end
  valid = any(T.rampMS == limits.rampMS & T.kernelCI > 0, 2);         % empty entries have zero for kernelCI
  if ~strcmp(limits.animal, 'All')
    valid = valid & sum(T.animal == limits.animal, 2);
  end
  if ~isempty(limits.oneDay)
    valid = valid & T.date == limits.oneDay;
  end
  valid = valid & (T.numStim > 0);                                    % only looking at sessions with stimulation
  valid = valid & (T.stimCorrects > limits.minTrials);                % enough correct trials
  valid = valid & T.stimFails > limits.minTrials;                    	% enough fail trials
  valid = valid & T.kernelPeak > limits.criterion;
  
  % d', no stim, and decrement limits
  if limits.minDPrime ~= -1
    valid = valid & double(T.noStimDPrime) >= limits.minDPrime;
  end
  if limits.minDec ~= -1
    valid = valid & double(T.noStimDPrime) - double(T.stimDPrime) >= limits.minDec;     % min behavioral decrement in hit rate
  end
  if limits.maxMeanPowerMW ~= -1
    valid = valid & T.meanPowerMW <= limits.maxMeanPowerMW;
  end
  U = T(valid, :);
%  	U.dPrime(U.dPrime == Inf) = NaN; 
	U.noStimDPrime(U.noStimDPrime == Inf) = NaN; 
	U.stimDPrime(U.stimDPrime == Inf) = NaN; 
 
  % Some limits related to over-session performance.  We can requie a minimum number of sessions for each animal/ramp,
  % and a minimum average delta-d'.
  if isempty(limits.oneDay)
    animals = unique(U.animal);
    rampDurs = unique(limits.rampMS);
    for r = 1:length(rampDurs)
      for a = 1:length(animals)
        animalRows = U.animal == animals{a} & U.rampMS == rampDurs(r);
        if limits.minSessions > 0 && sum(animalRows) < limits.minSessions
          U = U(~animalRows, :);
          continue;
        end
        if limits.minAvgDeltaDPrime >= 0 && nanmean(U.noStimDPrime(animalRows) - U.stimDPrime(animalRows)) < limits.minAvgDeltaDPrime
          U = U(~animalRows, :);
        end
      end
    end
  end
  if size(U, 1) == 0
    if length(limits.rampMS) == 1
      fprintf('getSubset: No valid sessions found for %d ms ramps, minTrials %d and decrement %.2f\n', ...
        limits.rampMS, limits.minTrials, limits.minDec);
    else
      fprintf('getSubset: No valid sessions found for multiple ramps, minTrials %d and decrement %.2f\n', ...
        limits.minTrials, limits.minDec);
    end
  end
end

function limits = setLimits(theSubset)
  % Set up the default selection criteria
  limits.minSessions = 0;                	% require at least n sessions for each animal
  limits.minTrials = 0;
  limits.criterion = 0;
  limits.minDec = -1;
  limits.minDPrime = -1;
  limits.minAvgDeltaDPrime = 0.10;
  limits.maxMeanPowerMW = 0.25;
  limits.animal = {'All'};
  limits.oneDay = [];
  switch theSubset
    case {'Unfiltered', 'unfiltered'}
      limits.minSessions = 0;
      limits.minDec = -1;
      limits.minDPrime = -1;
      limits.minAvgDeltaDPrime = -1;
      limits.maxMeanPowerMW = 0.25;
      limits.rampMS = [0, 500];
    case {'All', 'all', 'all steps', 'All Steps', 'All steps', 'all ramps', 'All ramps', 'All Ramps'}
      switch theSubset
        case {'All', 'all'}
          limits.rampMS = [0, 500];
        case {'all steps', 'All Steps', 'All steps'}
          limits.rampMS = 0;
        case {'all ramps', 'All ramps', 'All Ramps'}
          limits.rampMS = 500;
      end
    case {'Example', 'example'}
      limits.rampMS = 0;
      limits.animal = {'902'};
      limits.oneDay = '2019-10-10';
    case {'oneOff', 'oneoff', 'One Off'}
      limits.rampMS = 0;
      limits.animal = {'1462', '1463'};
    otherwise
      fprintf('getSessionTable: unrecognized table type ''%s''\n', theSubset);
      limits = [];
      return;
  end
end
