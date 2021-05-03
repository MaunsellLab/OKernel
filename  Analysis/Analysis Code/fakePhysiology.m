function fakePhysiology

    animalName = '1180';
    outFileName = 'fakePhysiology.mat';
    numTrials = 500;
    numChannels = 32;
    triggerChannel = 129;
    spikesPerMS = rand(1, numChannels) * 10.0 / 1000.0;                 % firing rates (per MS) for each channel
    trialDurS = 2.5;
    optoDurS = trialDurS + 1;
    pulseDurS = 0.025;
    zeroOptoS = 10.0;
    
    spikes = zeros(ceil(mean(spikesPerMS * 1000 * trialDurS)) * numChannels * numTrials * 2, 3);
    numSpikes = 0;
    
    for t = 1:numTrials
        trial.trialEnd.data = 0;
        trial.powerStimUW.data = 500;
        trial.trialStart.timeS = t * trialDurS;
        trial.trialEnd.timeS = trial.trialStart.timeS + trialDurS - 0.001;
        trial.stimulusOn.time = trial.trialEnd.timeS - trialDurS / 4.0;
        trial.zeroOptoFromPhotoDiode = floor(trial.stimulusOn.time - trial.trialStart.timeS) * 1000;
        
        numSpikes = numSpikes + 1;
        spikes(numSpikes, :) = [triggerChannel, 0, trial.trialStart.timeS];
        bins = trialDurS * 1000;
        random = rand(numChannels, bins);
        for c = 1:numChannels
            for b = 1:bins
                if random(c, b) < spikesPerMS(c)
                    numSpikes = numSpikes + 1;
                    spikes(numSpikes, :) = [c, 0, trial.trialStart.timeS + b * 0.001];
                end 
            end
        end
        optoBins = optoDurS / pulseDurS;
        trial.optoStepTimesMS = ((0:optoBins - 1) * pulseDurS - zeroOptoS) * 1000.0;
        powerMW = trial.powerStimUW.data * 2.0;
        trial.optoStepPowerMW = randi([0 1], 1, optoBins) * powerMW;
        
        % Here we insert extra spikes after a specific pattern of opto stimulation
        
        pattern = [powerMW, powerMW, 0, 0];                             % target stimulation pattern
        powers = trial.optoStepPowerMW;                                 % pointer to powers
        for b = 1:length(powers) - length(pattern) - 1                  % scan the entire stimulus
            patternMatch = pattern == powers(b:b + length(pattern) - 1);
            if sum(patternMatch) == length(pattern)
%                 spikeTimeS = trial.stimulusOn.time - trial.zeroOptoFromPhotoDiode / 1000 + (b + length(pattern)) * pulseDurS;
                spikeTimeS = trial.trialStart.timeS + (b - 1 + length(pattern)) * pulseDurS;
                for c = 1:numChannels                                   % add pattern driven spikes
                    numSpikes = numSpikes + 1;
                    spikes(numSpikes, :) = [c, 0, spikeTimeS];
                end
            end
        end

        trials(t) = trial; %#ok<AGROW>
    end
    
    % save the results
    taskEvents = {trials};
    taskSpikes = {spikes(1:numSpikes, :)};
    taskNames = {'OK'};
    save(['~/Desktop/', animalName, '/', outFileName], 'taskNames', 'taskEvents', 'taskSpikes');
end
