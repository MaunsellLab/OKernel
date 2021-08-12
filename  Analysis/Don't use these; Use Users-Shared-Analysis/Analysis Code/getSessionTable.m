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
  controls = controlSessions(T);                                  % logic array flagging control sesions
  valid = any(T.rampMS == limits.rampMS & T.kernelCI > 0, 2) & ~controls; % empty entries have zero for kernelCI
  if ~strcmp(limits.animal, 'All')
    valid = valid & T.animal == limits.animal;
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
 	U.dPrime(U.dPrime == Inf) = NaN; 
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
      fprintf('getSessionTable: No valid sessions found for ''%s'' for animal ''%s'' on %d ms ramps, minTrials %d and decrement %.2f\n', ...
        mode, limits.animal, limits.rampMS, limits.minTrials, limits.minDec);
    else
      fprintf('getSessionTable: No valid sessions found for ''%s'' for animal ''%s'' on multiple ramps, minTrials %d and decrement %.2f\n', ...
        mode, limits.animal, limits.minTrials, limits.minDec);
    end
  end
end

%%
function controls = controlSessions(T)

% control measurement sesssion.  These are contralateral stimulation
% controls that should be excluded from general analysis

  cSessions = {
    {'1218', {'2020-06-06', '2020-06-07', '2020-06-08', '2020-06-09', '2020-06-10', '2020-06-11', '2020-06-12', ...
              '2020-06-13', '2020-06-14', '2020-06-15', '2020-06-16'}},...
    {'1220', {'2020-06-22', '2020-06-23', '2020-06-24', '2020-06-25', '2020-06-26', '2020-06-27', '2020-06-28'}},...
    {'1257', {'2020-05-30', '2020-05-31', '2020-06-01', '2020-06-02', '2020-06-03', '2020-06-04', '2020-06-05', ...
              '2020-06-06', '2020-06-07', '2020-06-08', '2020-06-09', '2020-06-10', '2020-06-11', '2020-06-12', ...
              '2020-06-13', '2020-06-14'}},...
  };
  
  controls = false(height(T), 1);
  for a = 1:length(cSessions)
    for d = 1:length(cSessions{a}{2})
      controls = controls | (T.animal == cSessions{a}{1} & T.date == cSessions{a}{2}{d});
    end
  end
end

function controls = prePostControlSessions(T)

% control measurement sesssion.  These are a subset of normal stimulation
% sessions that were done immediately before and after the control stimulation
% session.  They are matched in number to the control stimulation sessions
% to keep the S/N balanced between the two

  cSessions = {
    {'1218', {'2020-06-01', '2020-06-02', '2020-06-03', '2020-06-04', '2020-06-05', ...
              '2020-06-21', '2020-06-22', '2020-06-23', '2020-06-24', '2020-06-25'}},...
    {'1220', {'2020-06-17', '2020-06-18', '2020-06-19', '2020-06-20', '2020-06-21'}},...
    {'1257', {'2020-05-25', '2020-05-26', '2020-05-27', '2020-05-28', '2020-05-29', ...
              '2020-06-21', '2020-06-22', '2020-06-23', '2020-06-24', '2020-06-25'}},...
  };
  
  controls = false(height(T), 1);
  for a = 1:length(cSessions)
    for d = 1:length(cSessions{a}{2})
      controls = controls | (T.animal == cSessions{a}{1} & T.date == cSessions{a}{2}{d});
    end
  end
end

function controls = controlTestSessions(T)

% control measurement sesssion. These are contralateral stimulation
% controls that should be excluded from general analysis.  This list is
% balanced to have the same number of sessions as the pre/post control
% measurements, to keep the S/N balanced between the two.

  cSessions = {
    {'1218', {'2020-06-07', '2020-06-08', '2020-06-09', '2020-06-10', '2020-06-11', ...
              '2020-06-12', '2020-06-13', '2020-06-14', '2020-06-15', '2020-06-16'}},...
    {'1220', {'2020-06-22', '2020-06-23', '2020-06-24', '2020-06-25', '2020-06-28'}},...
    {'1257', {'2020-06-04', '2020-06-05', '2020-06-06', '2020-06-07', '2020-06-09', ...
              '2020-06-10', '2020-06-11', '2020-06-12', '2020-06-13', '2020-06-14'}},...
  };
  
  controls = false(height(T), 1);
  for a = 1:length(cSessions)
    for d = 1:length(cSessions{a}{2})
      controls = controls | (T.animal == cSessions{a}{1} & T.date == cSessions{a}{2}{d});
    end
  end
end
