function trials = validateTrials(row, trials)
  %% validateTrials -- make sure that every trial has an EOT and meanPower so the indices will align
 
  trials = validateOnField(row, trials, 'trialEnd', -1);
  trials = validateOnField(row, trials, 'meanPowerMW', -1);
  trials = validateOnField(row, trials, 'trial', trials(1).trial);

  % it's a little harder to do the check on fields that are further nested
  numTrials = length(trials);
	trialStructs = [trials(:).trial];
	if isfield(trialStructs, 'pulseContrast')                               % some files have pulseContrast as well
    numPulseContrast = length([trialStructs(:).pulseContrast]);           % get pulse contrast applied on each trial
    if numPulseContrast < numTrials
%       fprintf('  %s %s: adding %d missing pulseContrasts\n', row.animal, row.date, numTrials - numPulseContrast);
      for t = 1:numTrials
        if isempty(trials(t).trial.pulseContrast)
%           fprintf('     trial %d of %d\n', t, numTrials);
          trials(t).trial.pulseContrast = -1;
        end
      end
    end  
  end
end

function trials = validateOnField(row, trials, fieldName, defaultValue)

  if isfield(trials, fieldName)
    numField = length([trials(:).(fieldName)]);                           % get reaction times on each trial
  else
    numField = 0;
  end
  numTrials = length(trials);
  if numField < numTrials
%     fprintf('  %s %s: adding %d missing %s fields\n', row.animal, row.date, numTrials - numField, fieldName);
    for t = 1:numTrials
    	if ~isfield(trials, fieldName) || isempty(trials(t).(fieldName))
%         fprintf('     trial %d of %d\n', t, numTrials);
        trials(t).(fieldName) = defaultValue;
      end
    end
  end
end
