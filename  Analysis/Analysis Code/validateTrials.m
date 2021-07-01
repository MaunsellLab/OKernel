function trials = validateTrials(trials)
  %% validateTrials -- make sure that every trial has an EOT and meanPower so the indices will align
 
  trials = atLeastOnePerTrial(trials, 'trialEnd', -1);
  trials = atLeastOnePerTrial(trials, 'meanPowerMW', -1);
  trials = atLeastOnePerTrial(trials, 'trial', trials(1).trial);

  % it's a little harder to do the check on fields that are further nested
  numTrials = length(trials);
	trialStructs = [trials(:).trial];
	if isfield(trialStructs, 'pulseContrast')                               % some files have pulseContrast as well
    numPulseContrast = length([trialStructs(:).pulseContrast]);           % get pulse contrast applied on each trial
    if numPulseContrast < numTrials
      for t = 1:numTrials
        if isempty(trials(t).trial.pulseContrast)
          trials(t).trial.pulseContrast = -1;
        end
      end
    end  
  end
  
  % ensure that there each trial has no more than one reactTimeMS and one and only one trialEnd
  for t = 1:length(trials)                               
    if length(trials(t).reactTimeMS) > 1
        trials(t).reactTimeMS = trials(t).reactTimeMS(1);   
    end
    if ~isfield(trials(t), 'trialEnd') || isempty(trials(t).trialEnd)
        trials(t).trialEnd = -1;
    elseif length(trials(t).trialEnd) > 1
        trials(t).trialEnd = trials(t).trialEnd(1);   
    end
  end
end

function trials = atLeastOnePerTrial(trials, fieldName, defaultValue)

  if isfield(trials, fieldName)
    numField = length([trials(:).(fieldName)]);                           % get number of entries for each trial
  else
    numField = 0;
  end
  numTrials = length(trials);
  if numField < numTrials
    for t = 1:numTrials
    	if ~isfield(trials, fieldName) || isempty(trials(t).(fieldName))
        trials(t).(fieldName) = defaultValue;
      end
    end
  end
end
