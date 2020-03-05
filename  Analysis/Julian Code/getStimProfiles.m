function stimProfiles = getStimProfiles(trials, startTimeMS, endTimeMS, normalize)
% getStimProfiles()
% Return a millisecond-resolution copy of the optogenetic stimulus profile from each of a set of trials. Start and 
% end times are specified by arguments.  The profile is constructed using the optoStepTimeMS and optoStepPowerMW 
% arrays. If "normalize" is true, the power profiles will be normalized between zero and one.

    numTrials = size(trials, 2);
    if numTrials == 0
        stimProfiles = [];
        return
    end
    stimProfiles = zeros(numTrials, endTimeMS - startTimeMS);
    for t = 1:numTrials
        pulseDurMS = trials(t).trial.pulseDurMS;
        startBin = find(trials(t).optoStepTimesMS > startTimeMS, 1);                      % find the first bin
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
    end
    % We have to normalize each trial individually, because the user might change the mean power
    % during the day and move it above the max for some trials
    if normalize
        for t = 1:numTrials
            meanValue = (max(stimProfiles(t, :)) + min(stimProfiles(t, :))) / 2.0;
            stimProfiles(t, stimProfiles(t, :) < meanValue) = 0;
            stimProfiles(t, stimProfiles(t, :) >= meanValue) = 1;
        end
    end
end
