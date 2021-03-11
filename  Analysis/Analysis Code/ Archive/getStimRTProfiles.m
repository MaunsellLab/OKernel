function stimProfiles = getStimRTProfiles(trials, startTimeMS, endTimeMS, normalize)
% getStimProfiles()
% Return a millisecond-resolution copy of the optogenetic stimulus profile from each of a set of trials. This function
% differ from getStimProfiles() in aligning the profiles to the reaction time on each trial (getStimProfiles aligns
% to stimulus onset. Start and end times relative to RT are specified by arguments. The profile is constructed using 
% the optoStepTimeMS and optoStepPowerMW arrays. If "normalize" is true, the power profiles will be normalized between
% zero and one.

    numTrials = size(trials, 2);
    if numTrials == 0
        stimProfiles = [];
        return
    end
    stimProfiles = zeros(numTrials, endTimeMS - startTimeMS);
    for t = 1:numTrials                                                                 % for each trial
        meanPower = trials(t).meanPowerMW;
        if meanPower == 0
            if normalize
                meanPower = 0.5;
            end 
            stimProfiles(t, :) = meanPower * ones(1, endTimeMS - startTimeMS);
            continue;
        end       
        pulseDurMS = trials(t).trial.pulseDurMS;                                        % duration of stim pulses
        
        trials(t).optoStepTimesMS = trials(t).optoStepTimesMS - trials(t).reactTimeMS;  % offset times to RT
 %{
The stimProfile is an integral number of pulses long.  We need to pad it out at the head and tail to make it 
have a length of (endTimeMS - startTimeMS).  This is all pretty simple to implement if the startBin isn't the
first bin.  That shouldn't really happen, but when it does, we need to take special measures.  This is 
because the first bin (uniquely) is not guaranteed to be pulseDurMS from the second bin (phase is randomized
on each trial by inserting a partial bin at the start of the sequence.  We simply start with the second bin, and
pad backward, assuming that the earlier values have the same value as the first bin
%}       
        startBin = find(trials(t).optoStepTimesMS > startTimeMS, 1);                  	% find the first bin
        startBin = max(startBin, 2);
        endBin = find(trials(t).optoStepTimesMS <= endTimeMS - pulseDurMS, 1, 'last');	% find last bin
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
%       Power might be set up and down on a per-trial basis(!), so we have to normalize each trial individually
        if normalize
            stimProfiles(t, stimProfiles(t, :) < meanPower) = 0;
            stimProfiles(t, stimProfiles(t, :) >= meanPower) = 1;
        end
    end
end
