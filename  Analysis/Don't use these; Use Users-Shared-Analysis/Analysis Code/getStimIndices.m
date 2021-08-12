function [stimIndices, trials] = getStimIndices(trials)
% Find the trials between the first and last stimulated trials, counting only those with full pulse contrasts
  meanPower = [trials(:).meanPowerMW];                        % get power applied for each trial                        
  if sum(meanPower) == 0                                      % no opto stimulation in this session
      stimIndices = [];
      trials = [];
      return;
  end
  trialStructs = [trials(:).trial];
  if isfield(trialStructs, 'pulseContrast')                   % get rid of any trials with reduced opto power
      stimIndices = meanPower > 0 & [trialStructs.pulseContrast] == 1;	% only trials with contrast == 1
  else
      stimIndices = meanPower > 0;
  end
  if sum(stimIndices) == 0                                    % no trials with full opto contrast
    trials = [];
    return;
  end
  % only consider the range between the first and last stimulated trials -- trim the vectors
  firstStimIndex = find(stimIndices > 0, 1);                 	% first stimulated trial
	lastStimIndex = find(stimIndices > 0, 1, 'last');          	% last stimulated trial
  stimIndices = stimIndices(firstStimIndex:lastStimIndex);
  trials = trials(firstStimIndex:lastStimIndex);
end