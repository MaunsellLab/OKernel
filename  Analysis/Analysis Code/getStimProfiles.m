function stimProfiles = getStimProfiles(trials, startTimeMS, endTimeMS, normalize, alignRT)
% getStimProfiles()
% Return a millisecond-resolution copy of the optogenetic stimulus profile from each of a set of trials. Start and 
% end times are specified by input arguments.  The profile is constructed using the optoStepTimeMS and optoStepPowerMW 
% arrays. If "normalize" is true, the power profiles will be normalized between zero and one. "alignRT" true aligns
% the profiles on the reaction time, rather than the stimulus on time.

% When we are RT aligned, we must check that we have enough time on each trial to span the full kernel.  This isn't
% a problem for proper RT alignment. For for early alignment, the RT could approach the start of the trial.
% The "prepend" code will fill in the blank with the value in the first bin, but this means that a few trials that
% do this filling on most of the kernel can dominate, causing an artificial, progressing offset.  To avoid this,
% we screen for these trials and eliminated them.
  
  if nargin < 5
    alignRT = false;
  end
  if alignRT
    numTrials = size(trials, 2);
    enoughBins = true(1, numTrials);
    for t = 1:numTrials  
      if trials(t).optoStepTimesMS(1) - trials(t).reactTimeMS > startTimeMS
        enoughBins(t) = false;
      end
    end
    trials = trials(enoughBins);
  end
  numTrials = size(trials, 2);
  if numTrials == 0
      stimProfiles = [];
      return
  end
  
  stimProfiles = zeros(numTrials, endTimeMS - startTimeMS);
  for t = 1:numTrials
    meanPower = trials(t).meanPowerMW;
    if meanPower == 0 && normalize
      stimProfiles(t, :) = 0.5 * ones(1, endTimeMS - startTimeMS);
      continue;
    end
    pulseDurMS = trials(t).trial.pulseDurMS;   
    if alignRT
        trials(t).optoStepTimesMS = trials(t).optoStepTimesMS - trials(t).reactTimeMS;  % offset times to RT
    end
%{
    The stimProfile is an integral number of pulses long.  We need to pad it out at the head and tail to make it 
    have a length of (endTimeMS - startTimeMS).  This is all pretty simple to implement if the startBin isn't the
    first bin.  That shouldn't really happen, but it does, and we need to take special measures.  This is mostly
    because the first bin (uniquely) is not guaranteed to be pulseDurMS from the second bin (phase is randomized
    on each trial by inserting a partial bin at the start of the sequence.
%}
    startBin = find(trials(t).optoStepTimesMS > startTimeMS, 1);                      % find the first bin
    startBin = max(startBin, 2);
    endBin = find(trials(t).optoStepTimesMS <= endTimeMS - pulseDurMS, 1, 'last');    % find last bin      
    stimProfile = repelem(trials(t).optoStepPowerMW(startBin:endBin), pulseDurMS);

    prepend = [];
    prepend(1:trials(t).optoStepTimesMS(startBin) - startTimeMS) = trials(t).optoStepPowerMW(startBin - 1);
    numAppendBins = endTimeMS - trials(t).optoStepTimesMS(endBin) - pulseDurMS;
    if numAppendBins > 0
      append = [];
      append(1:numAppendBins) = trials(t).optoStepPowerMW(endBin + 1);
      stimProfiles(t, :) = [prepend, stimProfile, append];
    else
        stimProfiles(t, :) = [prepend, stimProfile];
    end
    % We must normalize each trial individually, because the user might change the mean power
    % during the day and move it above the anticipated max for some trials. We divide mean power by
    % 50 because there might be a ramp at the start of the trial and we count those as full power so the
    % kernel doesn't get pulled to zero.
    
    if normalize
      stimProfiles(t, stimProfiles(t, :) < meanPower / 50) = 0;
      stimProfiles(t, stimProfiles(t, :) >= meanPower / 50) = 1;
    end
  end
end
