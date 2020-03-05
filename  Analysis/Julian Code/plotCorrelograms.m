function plotCorrelograms(plotIndex, kernel, plotStartMS, plotEndMS, dataFiles)

    numBins = plotEndMS - plotStartMS;
    numFiles = length(dataFiles);
    numTrials = 0;
    for f = 1:numFiles
        load(dataFiles{f});                                         % load the session
        meanPower = [trials(:).meanPowerMW]; %#ok<USENS>
        stimTrials = meanPower > 0;
        trialEnds = [trials(:).trialEnd];
        stimIndices = stimTrials & (trialEnds == 0 | trialEnds == 2);
        numTrials = numTrials + sum(stimIndices);
    end
    y = zeros(1, numBins * numTrials);
    index = 1;
    for f = 1:numFiles
        load(dataFiles{f});                                         % load the session
        meanPower = [trials(:).meanPowerMW];
        stimTrials = meanPower > 0;
        trialEnds = [trials(:).trialEnd];
        stimIndices = stimTrials & (trialEnds == 0 | trialEnds == 2);
        profiles = getStimProfiles(trials(stimIndices), plotStartMS, plotEndMS, false);
        for p = 1:size(profiles, 1)
            y(index:index + numBins - 1) = profiles(p, :);
            index = index + numBins;
        end
    end
    y(y < 0.01) = 0;
    y(y >= 0.01) = 1;
    
    maxLags = 75;
    [cStim, lags] = xcorr(y, maxLags);
    cStim = (cStim - mean(cStim(floor(maxLags + 26):end)));
    cStim = cStim / max(cStim);
    ax = subplot(4, 3, plotIndex + 1);
    plot(lags(maxLags:end), cStim(maxLags:end), 'r');
    hold on;
   	[cKernel, lags] = xcorr(kernel, 'normalized');
    plot(lags(numBins:numBins + maxLags), cKernel(numBins:numBins + maxLags), 'b');
    ax.Title.String = 'Kernel Autocorrel.';
    ax.XGrid = 'on';
    ax.YGrid = 'on';

    ax = subplot(4, 3, plotIndex);						% default axes are 0 to 1
    plot(lags(numBins:end), cKernel(numBins:end), 'b');    
    ax.Title.String = 'Kernel Autocorrel.';
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    maxLags = floor(length(cKernel) / 2);
    [cStim, lags] = xcorr(y, maxLags);
    cStim = (cStim - mean(cStim(floor(maxLags + 26):end)));
    cStim = cStim / max(cStim);
    hold on;
    plot(lags(maxLags:end), cStim(maxLags:end), 'r');
end