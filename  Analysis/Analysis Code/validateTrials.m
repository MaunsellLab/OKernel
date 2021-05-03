function trials = validateTrials(trials)
  %% validateTrials -- make sure that every trial has an EOT and meanPower so the indices will align
 
  trials = validateOneField(trials, 'trialEnd', -1);
  trials = validateOneField(trials, 'meanPowerMW', -1);
  trials = validateOneField(trials, 'trial', trials(1).trial);

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
end

function trials = validateOneField(trials, fieldName, defaultValue)

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
