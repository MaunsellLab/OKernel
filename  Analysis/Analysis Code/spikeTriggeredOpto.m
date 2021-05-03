function spikeTriggeredOpto

    animalName = '1180';
    fullFileName = '200320R2_MRG.mat';
%     fullFileName = 'fakePhysiology.mat';
	[~, name, ~] = fileparts(fullFileName);
    displayName = strrep(name, '_', '-');
    numChannels = 32;
    minSpikes = 100;
    minCIFactor = 4.0;
    startKernelMS = -350;
    endKernelMS = 100;
    triggerChannel = 129;

    % load and extract the data
    load(['~/Desktop/', animalName, '/', fullFileName], 'taskNames', 'taskEvents', 'taskSpikes');
    assert(strcmp(taskNames{1}, 'OK'), 'Data from file are not identified as coming from OKernel');    
    trials = taskEvents{1};
    spikes = taskSpikes{1};
    clear taskNames taskEvents taskSpikes
    
%     for c = 1:numChannels
    for c = 31:31
        channelIndices = spikes(:, 1) == c;
        spikeNum = sum(channelIndices);
        if spikeNum < minSpikes
            continue;
        end
        spikeCount = 0;
        triggerIndices = spikes(:, 1) == triggerChannel;
        profiles = zeros(spikeNum, endKernelMS - startKernelMS);
        for t = 1:length(trials)
            % only use correctly completed, opto stimulated trials
            if trials(t).trialEnd.data ~= 0 || trials(t).powerStimUW.data == 0
                continue;
            end

            % find the targeted channel spikes for this trial
            trialStartS = trials(t).trialStart.timeS;
            trialEndS = trials(t).trialEnd.timeS;
            validSpikes = channelIndices & spikes(:, 3) >= trialStartS & spikes(:, 3) <= trialEndS;
            if sum(validSpikes) == 0                            % no spikes, continue
                continue;
            end
            triggerEvents = triggerIndices & spikes(:, 3) >= trialStartS & spikes(:, 3) <= trialEndS;
            if sum(triggerEvents) ~= 1                          % one and only one photodiode trigger in trial
                fprintf('Warning: %d photodiode triggers on trial %d\n', sum(triggerEvents), t);
                continue;
            end
            optoStartS = spikes(triggerEvents, 3);
            
            % build the normalized optogenetic stimulus description for this trial in spike time base
            optoStepTimesS = (trials(t).optoStepTimesMS - trials(t).optoStepTimesMS(1)) / 1000.0 + optoStartS;
            assert(optoStepTimesS(end) > trialEndS, 'Trial extends beyond the opto stimulus');
            pulseDurMS = trials(t).optoStepTimesMS(3) - trials(t).optoStepTimesMS(2);
            % create a millisecond array of optogenetic power starting from the beginning of the opto stimulus train
            % the first bin is typically less than pulseDurMS long, because that is how we randomized opto pulse phase
            firstBinDurMS = trials(t).optoStepTimesMS(2) - trials(t).optoStepTimesMS(1);
            stimProfile = [];
            stimProfile(1:firstBinDurMS) = trials(t).optoStepPowerMW(1);
            endBin = length(trials(t).optoStepTimesMS);         % the last bin      
            stimProfile = [stimProfile repelem(trials(t).optoStepPowerMW(2:endBin), pulseDurMS)]; %#ok<AGROW>
            % We have to normalize each trial individually, because the user might change the mean power
            % during the day and move it above the max for some trials
            stimProfile(stimProfile(:) < trials(t).powerStimUW.data / 1000.0) = 0;
            stimProfile(stimProfile(:) >= trials(t).powerStimUW.data / 1000.0) = 1;

            % sum into the kernel using each spike time
            spikeTimesS = spikes(validSpikes, 3);
            for s = 1:length(spikeTimesS)
%                 firstStimBin = floor(((spikeTimesS(s) + startKernelMS / 1000.0) - optoStartS) * 1000.0);
                firstStimBin = ceil(((spikeTimesS(s) + startKernelMS / 1000.0) - optoStartS) * 1000.0);
                lastStimBin = firstStimBin + endKernelMS - startKernelMS - 1;
                % Check we have enough kernel.  We elimiante the first 200 ms because the power ramped up there
                if firstStimBin > 200 && lastStimBin <= length(stimProfile)
                    spikeCount = spikeCount + 1;
                    profiles(spikeCount, :) = stimProfile(firstStimBin:lastStimBin);
                    
%                     for i = 297:301
%                         fprintf(' %3d', stimProfile(firstStimBin + i));
%                     end
%                     fprintf('   ');
%                     for i = 322:326
%                         fprintf(' %3d', stimProfile(firstStimBin + i));
%                     end
%                     fprintf('   ');
%                     for i = 347:352
%                         fprintf(' %3d', stimProfile(firstStimBin + i));
%                     end
%                     fprintf(' %.3f\n', spikeTimesS(s) - optoStepTimesS(2));
                    
                    
                end
            end
        end

        % only display plots that are statistically reliable
        CI = stimCI(spikeCount);
        averageOpto = mean(profiles(1:spikeCount, :));
        if min(averageOpto) > 0.5 - minCIFactor * CI && max(averageOpto) < 0.5 + minCIFactor * CI
            continue;
        end
        
        %plot it
        axisHandle = figure(1);        
        set(axisHandle, 'Units', 'inches', 'Position', [25.5, 7.5, 8.5, 11]);
        clf;
        plot(averageOpto);
        xlim([1, size(profiles, 2)]);
        hold on;
        h = fill([0, endKernelMS - startKernelMS, endKernelMS - startKernelMS, 0], ...
            [0.5 + CI, 0.5 + CI, 0.5 - CI, 0.5 - CI], [0.8, 0.8, 0.8]);
        set(h, 'linestyle', ':', 'facealpha', 0.25);
        a = axis();
        plot([-startKernelMS, -startKernelMS], [a(3), a(4)], 'k:');
        plot([a(1), a(2)], [0.5, 0.5], ':', 'color', [0.0, 0.5, 1.0]);
        tickValues = fix(startKernelMS / 100) * 100 : 100 : fix(endKernelMS / 100) * 100;
        tickLabels = cell(1, length(tickValues));
        for t = 1:length(tickValues)
            tickLabels{t} = sprintf('%d', tickValues(t));
        end
        ax = gca;
        ax.XTick = (tickValues - startKernelMS);
        ax.XTickLabel = tickLabels;
        ax.XGrid = 'on';
        ax.XAxis.MinorTickValues = tickValues - startKernelMS + 50;
        ax.XMinorGrid = 'on';

        headerText = cell(1, 1);
        headerText{1} = sprintf('Spike Triggered Average');
        headerText{length(headerText) + 1} = sprintf('Animal %s, File %s', animalName, displayName);
        headerText{length(headerText) + 1} = sprintf('Spike Channel %d ', c);
        headerText{length(headerText) + 1} = sprintf('%d Spikes, %d Trials', spikeCount, length(trials));
        text(a(1) + 0.05 * (a(2) - a(1)), a(4) - 0.05 * (a(4) - a(3)), headerText, 'VerticalAlignment', 'top');
        axis(a);
%         saveas(gcf,['~/Desktop/', animalName, '/', fullFileName, sprintf('-%d.pdf', c)]);
    end
end