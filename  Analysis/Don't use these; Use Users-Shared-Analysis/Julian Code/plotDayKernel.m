function plotDayKernel(~, file, trials)
    % if we've already plotted some trials and this one isn't stimulated or usuable, do nothing;
    if figureHasSubplotTitleString(1, 'Kernel')         % return if we've already
        if (trials(end).trialEnd ~= 0 && trials(end).trialEnd ~= 2) || trials(end).meanPowerMW == 0
            return
        end
    end
    startTimeMS = -200;
    endTimeMS = 200;
    meanPowers = [trials(:).meanPowerMW];
    stimIndices = meanPowers > 0;
    trialEnds = [trials(:).trialEnd];
    trialProfile = getStimProfiles(trials(end), startTimeMS, endTimeMS, false);
    doOneKernelPlot(9, trialProfile, startTimeMS, endTimeMS, 'Optical Stim Last Trial', 'Power (mW)', 0);
    hitIndices = stimIndices & trialEnds == 0;
    missIndices = stimIndices & trialEnds == 2;
    % If we have kernelRTMax and Min, move excluded RTs from hits to misses
    if isfield(file, 'kernelRTMinMS') && isfield(file, 'kernelRTMaxMS')
        RTs = [trials(:).reactTimeMS];
        earlyIndices = hitIndices & RTs < file.kernelRTMinMS;
        lateIndices = hitIndices & RTs >= file.kernelRTMaxMS;
        hitIndices = hitIndices & ~(earlyIndices | lateIndices);
        missIndices = missIndices | lateIndices;
    end
    % Use the hit and miss indices to construct the profiles
    if sum(hitIndices) > 0
        hitProfiles = getStimProfiles(trials(hitIndices), startTimeMS, endTimeMS, true);
        hitMean = mean(hitProfiles);
    else
        hitMean = [];
    end
    if sum(missIndices) > 0
        missProfiles = getStimProfiles(trials(missIndices), startTimeMS, endTimeMS, true);
        missMean = mean(missProfiles);
    else
        missMean = [];
    end
    % If this is a correct trial, or there is no correct plot yet, do that plot
    if trials(end).trialEnd == 0 || ~figureHasSubplotTitleString(1, 'Kernel Hit Trials')
        if ~isempty(hitMean)
            hitCI = stimCI(size(hitProfiles, 1));
            plotTitle = sprintf('Hit Kernel (n=%d)', sum(hitIndices));
            doOneKernelPlot(10, hitMean, startTimeMS, endTimeMS, plotTitle, 'Normalized Power', hitCI);
        end
    end
    % If this is a miss trial, or there is no miss plot yet, do that plot
    if trials(end).trialEnd == 2 || ~figureHasSubplotTitleString(1, 'Kernel Miss Trials')
        if ~isempty(missMean)
            missCI = stimCI(size(missProfiles, 1));
            plotTitle = sprintf('Miss Kernel (n=%d)', sum(missIndices));
            doOneKernelPlot(11, missMean, startTimeMS, endTimeMS, plotTitle, 'Normalized Power', missCI);
        end
    end
    if ~isempty(hitMean) && ~isempty(missMean)
        plotTitle = sprintf('Total Kernel (n=%d)', sum(hitIndices) + sum(missIndices));
        totalMean = hitMean - missMean;
        doOneKernelPlot(12, totalMean, startTimeMS, endTimeMS, plotTitle, 'Normalized Power', sqrt(hitCI^2 + missCI^2));
        % make sure that the hit and miss kernels (which both exist) are on the same y-axis
        subplot(4, 3, 10);
        ax = gca;
        y0 = ylim;
        subplot(4, 3, 11);
        y1 = ylim;
        theLimits = [min(y1(1), y0(1)), max(y1(2), y0(2))];
        ylim(gca, theLimits);
        ylim(ax, theLimits);
    end