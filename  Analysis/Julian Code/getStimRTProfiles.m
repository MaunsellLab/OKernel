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
        pulseDurMS = trials(t).trial.pulseDurMS;                                        % duration of stim pulses
        trials(t).optoStepTimesMS = trials(t).optoStepTimesMS - trials(t).reactTimeMS;  % offset times to RT
        startBin = find(trials(t).optoStepTimesMS > startTimeMS, 1);                	% find the first bin
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
    end
    if normalize
        meanValue = (max(max(stimProfiles)) + min(min(stimProfiles))) / 2.0;
        stimProfiles(stimProfiles < meanValue) = 0;
        stimProfiles(stimProfiles >= meanValue) = 1;
    end
end
