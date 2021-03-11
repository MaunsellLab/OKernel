function test

    load('/Users/maunsell/Desktop/1180/200320R2_MRG.mat', 'taskSpikes', 'taskEvents');
    spikes = taskSpikes{1};
    trials = taskEvents{1};
    indices = spikes(:,1) == 31;
    spikeTimes = spikes(indices, 3);
	triggerTimesS = spikes(spikes(:, 1) == 129, 3);
    numTrials = length(trials);
%     powerIndices = powerUW > 0;
    powerUW = zeros(1, numTrials);
	startTimeS = zeros(1, numTrials);
	endTimeS = zeros(1, numTrials);

    for t = 1:numTrials
        powerUW(t) = trials(t).powerStimUW.data;
        startTimeS(t) = trials(t).trialStart.timeS;
        endTimeS(t) = trials(t).trialEnd.timeS;
        trialSpikeTimes = spikeTimes(spikeTimes >= startTimeS(t) & spikeTimes <= endTimeS(t));
        if powerUW(t) == 0 || isempty(trialSpikeTimes)
            continue;
        end
        triggerTimeS = triggerTimesS(find(triggerTimesS > startTimeS(t), 1));
        if isempty(triggerTimeS) || triggerTimeS > endTimeS(t)
            continue;
        end
        formatStr = '%9.3f';
        profile = zeros(100000, 7);
        pCount = 0;
        optoTimeS = (trials(t).optoStepTimesMS - trials(t).optoStepTimesMS(1)) / 1000.0 + triggerTimeS;
    	fprintf('Trial %d, Trigger %.3f:\n', t, triggerTimeS);
        for i = 1:length(trialSpikeTimes)
            alignBin = find(optoTimeS < trialSpikeTimes(i), 1, 'last');
            pCount = pCount + 1;
            profile(pCount, :) = trials(t).optoStepPowerMW(alignBin - 3:alignBin + 3);
            fprintf('\nEvent:%43.3f\n', trialSpikeTimes(i));
            fprintf(' Times:      ');
            for b = alignBin - 3:alignBin + 3
                fprintf(formatStr, optoTimeS(b));
            end
            fprintf('\n');
            fprintf(' Diff:       ');
            for b = alignBin - 3:alignBin + 3
                fprintf(formatStr, optoTimeS(b) - trialSpikeTimes(i));
            end
            fprintf('\n');
            fprintf(' Value:      ');
            for b = alignBin - 3:alignBin + 3
                fprintf(formatStr, trials(t).optoStepPowerMW(b));
            end
            fprintf('\n');
        end
        mean(profile(1:pCount, :))
    end
%     spikesPerTrial = zeros(1, numTrials);
%     for t = 1:numTrials - 1
%         spikesPerTrial(t) = sum(spikeTimes >= startTimeS(t) & spikeTimes < startTimeS(t + 1));
%     end
%     spikesPerTrial(numTrials) = sum(spikeTimes > startTimeS(numTrials));
%     theTrial = 3;
% 	tIndices = spikes(:,1) == 129 & spikes(:,3) >= startTimeS(theTrial) & spikes(:,3) < endTimeS(theTrial);
%     fprintf('Trial %d: Start %.3f, End %.3f\n', theTrial, startTimeS(theTrial), endTimeS(theTrial));
%     fprintf('         Trigger %.3f\n', spikes(tIndices, 3));
%     fprintf('         Opto times   %6.3f %6.3f %6.3f %6.3f\n', trials(theTrial).optoStepTimesMS(1:4) / 1000.0);
%     fprintf('         Offset times %6.3f %6.3f %6.3f %6.3f\n', (trials(theTrial).optoStepTimesMS(1:4) - trials(theTrial).optoStepTimesMS(1))...
%         / 1000.0 + spikes(tIndices, 3));
%     
%     theSpikeTimes = spikeTimes(spikeTimes >= startTimeS(theTrial) & spikeTimes < endTimeS(theTrial));
%     fprintf('         Artifacts ');
%     for a = 1:length(theSpikeTimes)
%         fprintf(' %6.3f', theSpikeTimes(a));
%     end
%     fprintf('\n');
%     fprintf('         Framing: ');
%     for a = 1:length(theSpikeTimes)
%         fprintf(' %6.3f', mod(theSpikeTimes(a) - trials(theTrial).optoStepTimesMS(2) / 1000.0, 0.025));
%     end
%     fprintf('\n');

end
