function stimProfiles = getStimRTProfiles(trials, startTimeMS, endTimeMS)
% getStimRTProfiles()
% Return a millisecond-resolution copy of the optogenetic stimulus profile
% from each of a set of trials. This is a special case of getStimProfiles
% that stretches or compresses time so that the period between stimOn and
% RT is fixed across all trials
    
  numTrials = size(trials, 2);
  if numTrials == 0
    stimProfiles = [];
    return
  end
  
  [stimToRTMS, postRTMS] = stimRTLimits();
	preStimMS = -((endTimeMS - startTimeMS) - stimToRTMS - postRTMS);     % start relative to stim on
  kernelMS = endTimeMS - startTimeMS;                                   % duration of kernel
  stimProfiles = zeros(numTrials, kernelMS);
  
  for t = 1:numTrials
    meanPower = trials(t).meanPowerMW;
    if meanPower == 0
      meanPower = 0.5;
      stimProfiles(t, :) = meanPower * ones(1, kernelMS);
      continue;
    end
    pulseDurMS = trials(t).trial.pulseDurMS;
    
%     trials(t).optoStepTimesMS = trials(t).optoStepTimesMS - trials(t).reactTimeMS;  % offset times to RT
%{
    The stimProfile is an integral number of pulses long.  We need to pad it out at the head and tail to make it 
    have a length of (endTimeMS - startTimeMS).  This is all pretty simple to implement if the startBin isn't the
    first bin.  That shouldn't really happen, but when it does, we need to take special measures.  This is mostly
    because the first bin (uniquely) is not guaranteed to be pulseDurMS from the second bin (phase is randomized
    on each trial by inserting a partial bin at the start of the sequence.
%}

    startBin = find(trials(t).optoStepTimesMS > preStimMS, 1);                          % find the first bin
    startBin = max(startBin, 2);
    % We need to move endBin so that it is at a time after the RT
%     endBin = find(trials(t).optoStepTimesMS <= endTimeMS - pulseDurMS, 1, 'last');    % find last bin 
    if trials(t).reactTimeMS + postRTMS > trials(t).optoStepTimesMS(end)
      fprintf('You''re looking for a time past the end of the opto times\n');
    end
    endBin = find(trials(t).optoStepTimesMS <= trials(t).reactTimeMS + postRTMS, 1, 'last');
    baseProfile = repelem(trials(t).optoStepPowerMW(startBin:endBin), pulseDurMS);
    
    prepend = [];
    prepend(1:trials(t).optoStepTimesMS(startBin) - preStimMS) = trials(t).optoStepPowerMW(startBin - 1);
    numAppendBins = trials(t).reactTimeMS + postRTMS - trials(t).optoStepTimesMS(endBin);
    if numAppendBins > 0
      append = [];
      append(1:numAppendBins) = trials(t).optoStepPowerMW(endBin + 1);
      rawProfile = [prepend, baseProfile, append];
    else
      rawProfile = [prepend, baseProfile];
    end
    
    % Use resample to get the profile to the correct length. Because the
    % profile is nothing but transients, this causes a bit of ringing, but
    % we will eliminate that when we normalize (below).
    stimProfile = resample(rawProfile, kernelMS, length(rawProfile), 5, 25);
    % We have to normalize each trial individually, because the user might change the mean power
    % during the day and move it above the max for some trials
    stimProfile(stimProfile(:) < meanPower) = 0;
    stimProfile(stimProfile(:) >= meanPower) = 1;
    stimProfiles(t, :) = stimProfile;
  end
end
