function aggregate_JDC()

% The plots should include:
%   1) Autocorrelogram of stimulus -- should be the triangle with 50 ms base (SEM for plots?)
%   2) Autocorrelogram of kernel -- should be broader than the stim auto-corr.
%   3) Cross-correlogram of stimulus-kernel (trace of N/P kernel?)

[IncludeDays] = makeIncludeDays;


combineAcrossAnimals = 1;
nowsubject = '1112';
stepStim = 0;
rampStim = 1;

%datapath for development computer
datapath = '/Users/julian/Documents/MATLAB/OKernel/';

%datapath for John's computer
% datapath = '/Users/Shared/Data/OKernel/';

dataFiles = {};
if stepStim
    if combineAcrossAnimals
        for subj = 2:size(IncludeDays,2)
            if ~isempty(IncludeDays{2,subj})
                for i = 1:size(IncludeDays{2,subj},1)
                dataFiles = {dataFiles{:},...
                    [datapath IncludeDays{1,subj} '/MatFiles/' IncludeDays{2,subj}(i,:) '.mat']};
                end
            end
            
            
        end
    else
        for subj = 2:size(IncludeDays,2)
            if strcmp(nowsubject,IncludeDays{1,subj})
                for i = 1:size(IncludeDays{2,subj},1)
                dataFiles = {dataFiles{:},...
                    [datapath IncludeDays{1,subj} '/MatFiles/' IncludeDays{2,subj}(i,:) '.mat']};
                end
            end
            
            
        end
    end
end

if rampStim
    if combineAcrossAnimals
        for subj = 2:size(IncludeDays,2)
            if ~isempty(IncludeDays{3,subj})
                for i = 1:size(IncludeDays{3,subj},1)
                dataFiles = {dataFiles{:},...
                    [datapath IncludeDays{1,subj} '/MatFiles/' IncludeDays{3,subj}(i,:) '.mat']};
                end
            end
            
            
        end
    else
        for subj = 2:size(IncludeDays,2)
            if strcmp(nowsubject,IncludeDays{1,subj})
                for i = 1:size(IncludeDays{3,subj},1)
                dataFiles = {dataFiles{:},...
                    [datapath IncludeDays{1,subj} '/MatFiles/' IncludeDays{3,subj}(i,:) '.mat']};
                end
            end
            
            
        end
    end
end

    
    numFiles = length(dataFiles);
    pathParts = split(dataFiles.', '/');
    if numFiles == 1
        numAnimals = 1;
    else
        animals = cell(1, numFiles);
        dates = cell(1, numFiles);
        for a = 1:numFiles
            animals{a} = pathParts{a, 6};
            dates{a} = pathParts{a, 8};
        end
        numAnimals = length(unique(animals));
    end
    plotStartMS = -400;
    plotEndMS = 400;
    plotRTStartMS = -600;
    plotRTEndMS = 200;
    bins = plotEndMS - plotStartMS;
    kernelRTMinMS = 250;
    kernelRTMaxMS = 700;
    % allocate buffers
    dayHitSums = zeros(numFiles, bins);
    dayNumHits = zeros(1, numFiles);
    dayHitScrambleSums = zeros(numFiles, bins);
    dayNumScrambleHits = zeros(1, numFiles);
    dayRTSums = zeros(numFiles, bins);
    dayNumRT = zeros(1, numFiles);
    dayMissSums = zeros(numFiles, bins);
    dayNumMisses = zeros(1, numFiles);
    dayMissScrambleSums = zeros(numFiles, bins);
    dayNumScrambleMisses = zeros(1, numFiles);
    dayFASums = zeros(numFiles, bins);
    dayNumFA = zeros(1, numFiles);    
    correctRTs = [];
    wrongRTs = [];
    missRTs = [];
    maxRespTimeMS = -1000;
    minRespTimeMS = 1000;
    % process each session (file)
    for f = 1:numFiles
        if numAnimals > 1
            fprintf('%s %s\n', animals{f}, dates{f});
        end
        load(dataFiles{f});                                         % load the session
        minRespTimeMS = min(minRespTimeMS, file.tooFastMS);         % set min/max response times
        maxRespTimeMS = max(maxRespTimeMS, file.rewardedLimitMS);   
        meanPower = [trials(:).meanPowerMW];                        % get power for each trial                        
        stimIndices = meanPower > 0;                                % get trials with opto stimulus
        trialEnds = [trials(:).trialEnd];                           % get eotCodes
        trialStructs = [trials(:).trial];                           % get trial structs
        preStimMS = [trialStructs(:).preStimMS];
        for t = 1:length(trials)                                    % allow only one RT per trial
            if length(trials(t).reactTimeMS) > 1
                trials(t).reactTimeMS = trials(t).reactTimeMS(1);   %#ok<AGROW>
            end
        end
        RTs = [trials(:).reactTimeMS];                              % get all RTs
        hitIndices = stimIndices & trialEnds == 0;                  % get different trial end types
        faIndices = (stimIndices & trialEnds == 1) & (preStimMS - 250 + plotStartMS > RTs);
        earlyIndices = hitIndices & RTs < kernelRTMinMS;            % hits outside of the allowed RT range
        lateIndices = hitIndices & RTs >= kernelRTMaxMS;
        hitIndices = hitIndices & ~(earlyIndices | lateIndices);    % strip out hits outside of range
        faIndices = faIndices | earlyIndices;                       % put early hits into FAs
        missIndices = (stimIndices & trialEnds == 2) | lateIndices;	% put late hits into misses
        % Use the hit and miss indices to construct the profiles
        numHits = sum(hitIndices);
        if numHits > 0
            profiles = getStimProfiles(trials(hitIndices), plotStartMS, plotEndMS, true);
            [dayNumHits(f), dayHitSums(f, :)] = addToSums(profiles);
            profiles = getStimRTProfiles(trials(hitIndices), plotRTStartMS, plotRTEndMS, true);
            [dayNumRT(f), dayRTSums(f, :)] = addToSums(profiles);
        end
        numMisses = sum(missIndices);
        if numMisses > 0
            profiles = getStimProfiles(trials(missIndices), plotStartMS, plotEndMS, true);
            [dayNumMisses(f), dayMissSums(f, :)] = addToSums(profiles);
        end
        if numHits + numMisses > 0
            validTrials = hitIndices | missIndices;
            randomList = randperm(numHits + numMisses);
            validIndices = find(validTrials > 0);
            scrambledHits = false(1, length(validTrials));
            scrambledMisses = false(1, length(validTrials));
            for i = 1:length(validIndices)
                if randomList(i) <= numHits
                    scrambledHits(validIndices(i)) = 1;
                else
                    scrambledMisses(validIndices(i)) = 1;
                end
            end
            profiles = getStimProfiles(trials(scrambledHits), plotStartMS, plotEndMS, true);
            [dayNumScrambleHits(f), dayHitScrambleSums(f, :)] = addToSums(profiles);
            profiles = getStimProfiles(trials(scrambledMisses), plotStartMS, plotEndMS, true);
            [dayNumScrambleMisses(f), dayMissScrambleSums(f, :)] = addToSums(profiles);
        end
        if sum(faIndices) > 0
            profiles = getStimProfiles(trials(faIndices), plotStartMS, plotEndMS, true);
            [dayNumFA(f), dayFASums(f, :)] = addToSums(profiles);
        end 
        
        % Add to the RT distributions
        correctIndices = find([trials(:).trialEnd] == 0);      % correct trials
        if ~isempty(correctIndices)
            correctRTs = [correctRTs [trials(correctIndices).reactTimeMS]]; %#ok<AGROW>
        end
        faIndices = find([trials(:).trialEnd] == 1);           % false alarm trials
        if ~isempty(faIndices)
            wrongRTs = [wrongRTs [trials(faIndices).reactTimeMS]]; %#ok<AGROW>
        end
        missIndices = find([trials(:).trialEnd] == 2);         % miss trials    
        if ~isempty(missIndices)
            missRTs = [missRTs [trials(missIndices).reactTimeMS]]; %#ok<AGROW>
            missRTs(missRTs < 0) = 100000;        % include misses in count, but don't let them display on plot
        end
    end
 
    % Get the average of session kernels
    hitMean = zeros(1, bins); 
    missMean = zeros(1, bins); 
    hitScrambleMean = zeros(1, bins); 
    missScrambleMean = zeros(1, bins); 
    faNormMean = zeros(1, bins); 
    rtMean = zeros(1, bins); 
    for f = 1:numFiles
        hitMean = hitMean + dayHitSums(f, :) / dayNumHits(f);  % sum in the daily averages
        missMean = missMean + dayMissSums(f, :) / dayNumMisses(f);
        hitScrambleMean = hitScrambleMean + dayHitScrambleSums(f, :) / dayNumScrambleHits(f);  % sum in the daily averages
        missScrambleMean = missScrambleMean + dayMissScrambleSums(f, :) / dayNumScrambleMisses(f);
        if dayNumFA(f) > 0
            faNormMean = faNormMean + dayFASums(f, :) / dayNumFA(f);
        end
        rtMean = rtMean + dayRTSums(f, :) / dayNumRT(f);
    end
    hitMean = hitMean / numFiles;
    missMean = missMean / numFiles;
    hitScrambleMean = hitScrambleMean / numFiles;
    missScrambleMean = missScrambleMean / numFiles;
    faNormMean = faNormMean / sum(dayNumFA > 0);
    rtMean = rtMean / numFiles;
    
    % Plot the results
   	h = figure(1);
    set(h, 'Units', 'inches', 'Position', [20, 5, 8.5, 11]);
    clf;
    axisHandle = subplot(4, 3, 1);						% default axes are 0 to 1
    set(axisHandle, 'Visible', 'off');
    set(axisHandle, 'OuterPosition', [0.02 0.75, 0.25, 0.2]);
    text(0.00, 1.25, 'OKernel', 'FontWeight', 'bold', 'FontSize', 16);
    headerText = cell(1, 1);
    headerText{1} = sprintf('%d sessions from %d animals', numFiles, numAnimals);
    headerText{length(headerText) + 1} = sprintf('Reaction times from %d to %d ms', kernelRTMinMS, kernelRTMaxMS);
    text(0.00, 1.00, headerText, 'VerticalAlignment', 'top');
    if ~isempty(hitMean)
        numHits = sum(dayNumHits);
        plotTitle = sprintf('Hit Kernel (n=%d)', numHits);
        hitCI = stimCI(numHits);
        doOneKernelPlot(hitMean, plotStartMS, plotEndMS, plotTitle, 4, 'Normalized Power', hitCI);
    end
    if ~isempty(missMean)
        numMisses = sum(dayNumMisses);
        plotTitle = sprintf('Miss Kernel (n=%d)', numMisses);
        missCI = stimCI(sum(dayNumMisses));
        doOneKernelPlot(missMean, plotStartMS, plotEndMS, plotTitle, 5, 'Normalized Power', missCI);
    end
    if ~isempty(hitMean) && ~isempty(missMean)
        plotTitle = sprintf('Total Kernel (n=%d)', numHits + numMisses);
        total = hitMean - missMean;
        doOneKernelPlot(total, plotStartMS, plotEndMS, plotTitle, 6, 'Normalized Power', sqrt(hitCI^2 + missCI^2));
    end
    if ~isempty(hitScrambleMean) && ~isempty(missScrambleMean)
        numScrambleHits = sum(dayNumScrambleHits);
        numScrambleMisses = sum(dayNumScrambleMisses);
        plotTitle = sprintf('Scrambled Kernel (n=%d)', numScrambleHits + numScrambleMisses);
      	total = hitScrambleMean - missScrambleMean;
        doOneKernelPlot(total, plotStartMS, plotEndMS, plotTitle, 9, 'Normalized Power', sqrt(hitCI^2 + missCI^2));
    end
    if ~isempty(rtMean)
        plotTitle = sprintf('RT Aligned (n=%d)', sum(dayNumRT));
        doOneKernelPlot(rtMean, plotRTStartMS, plotRTEndMS, plotTitle, 7, 'Normalized Power', stimCI(sum(dayNumRT)));
    end
    if ~isempty(faNormMean)
        plotTitle = sprintf('Norm. Avg. FAs (n=%d)', sum(dayNumFA));
        doOneKernelPlot(faNormMean, plotStartMS, plotEndMS, plotTitle, 8, 'Normalized Power', stimCI(sum(dayNumFA)));
    end
    sameYAxisScaling(4, 3, [4, 5, 7, 8]);
    sameYAxisScaling(4, 3, [6, 9]);
    RTHistogram(correctRTs, wrongRTs, missRTs, kernelRTMinMS, kernelRTMaxMS, minRespTimeMS, maxRespTimeMS);
    RTPDF(correctRTs, wrongRTs, missRTs, kernelRTMinMS, kernelRTMaxMS, minRespTimeMS, maxRespTimeMS);
%     plotCorrelograms(10, hitMean - missMean, plotStartMS, plotEndMS, dataFiles);
end

%%
function [numTrials, normSums] = addToSums(profiles)
%
    numTrials = size(profiles, 1);                                  % hits on this session
    meanPower = (max(max(profiles)) + min(min(profiles))) / 2.0;
    profiles(profiles < meanPower) = 0;
    profiles(profiles >= meanPower) = 1;
    normSums = sum(profiles);                                       % summed normalized profiles
end

% %%
% function powerFactor = checkPower(file, trials)
% % There were a few early files where the values saved were not the power, but the voltage that was
% % sent to the NIDAQ.  The difference was only about a factor of 1.5, but it needs to be corrected if
% % the plots are to be accurate in terms of power.
% %
% % NB: It later became clear that we needed to normalize the power of the stimuli across sessions to combine the data,
% % so this function got mooted.
% 
%     if isfield(file, 'text')
%         powerFactor = 1.0;
%     else
%         meanVs = [trials(:).meanPowerMW];               % it was called mW, but it was the volts out
%         stimIndex = find(meanVs > 0.0, 1);
%         maxV = max(trials(stimIndex).optoStepPowerMW);
%         maxPower = max([trials().meanPowerMW]) * 2.0;    % actual max power when stimulating
%         powerFactor = maxPower / maxV;
%     end
% end

%%
function doOneKernelPlot(profile, startTimeMS, endTimeMS, plotTitle, subIndex, label, theCI)

    subplot(4, 3, subIndex);
    cla;
    if contains(plotTitle, 'Total Kernel') || contains(plotTitle, 'Scrambled')
        theMean = 0;
    else
        theMean = 0.5;
    end
    bins = size(profile, 2);
    if nargin >= 7
        posCI = theMean + theCI;
        negCI = theMean - theCI;
        h = fill([0, bins, bins, 0], [posCI, posCI, negCI, negCI], [0.8, 0.8, 0.8]);
        set(h, 'linestyle', ':', 'facealpha', 0.25)
    end
    hold on;
    plot(profile, 'b');
    ax = gca;
    ax.XGrid = 'on';
    xlim(ax, [0, bins]);
    if contains(plotTitle, 'RT Aligned')
        set(gca,'XTick', [0, -startTimeMS, bins]);
        set(gca, 'XTickLabel', {sprintf('%d', startTimeMS), '0', sprintf('%d', endTimeMS)});
        xlabel('Time Relative to RT');
    else
        set(gca,'XTick', [0, -startTimeMS, -startTimeMS + 50 -startTimeMS + 100, bins]);
        set(gca, 'XTickLabel', {sprintf('%d', startTimeMS), '0', '', '100', sprintf('%d', endTimeMS)});
        xlabel('Time Relative to Stim On');
    end
    title(plotTitle);
    ylabel(label);
    plot([0 bins], [theMean, theMean], 'k:');
    hold off;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RT Histogram

 function RTHistogram(correctRTs, wrongRTs, missRTs, kernelRTMinMS, kernelRTMaxMS, minRespTimeMS, maxRespTimeMS)
 
    subplot(4, 3, 2);
    cla;
    hold on;
    timeLimit = 1000;
    RTBins = 100;
    while true
        edges = linspace(-1000, timeLimit, RTBins);
        nCorrect = histc(correctRTs, edges); %#ok<*HISTC>
        nWrong = histc(wrongRTs, edges);
        nMiss = histc(missRTs, edges);
        if max([nWrong, nCorrect, nMiss] <= 50)         % re-bin on next plot?
            break;
        else
            RTBins = min([RTBins * 2, 100]);
        end
    end
    if sum(nCorrect) + sum(nWrong) + sum(nMiss) > 0
        binSize = edges(2) - edges(1);
        bH = bar(edges + binSize / 2, [nCorrect(:), nMiss(:), nWrong(:)], 'stacked');
        set(bH, 'barWidth', 1, 'lineStyle', 'none');
        set(bH(1), 'FaceColor', [0 0 0.6]);
        set(bH(2), 'FaceColor', [0.6 0 0]);
        set(bH(3), 'faceColor', [0.6 0 0]);
        yLimits = get(gca, 'YLim');                % vertical line at stimulus on
        plot([0 0], yLimits, 'k');
        plot(minRespTimeMS * [1 1], yLimits, '--', 'Color', 0.5 * [0 1 0]);
        plot(maxRespTimeMS * [1 1], yLimits, '--', 'Color', 0.5 * [1 0 0]);
        plot(double(kernelRTMinMS) * [1 1], yLimits, ':', 'Color', 0.5 * [0 1 0]);
        plot(double(kernelRTMaxMS) * [1 1], yLimits, ':', 'Color', 0.5 * [1 0 0]);
    end
    set(gca, 'XLim', [-1000 timeLimit]);
    xlabel('Time Relative to Stimulus');
    title('Reaction Times');
end
 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RT PDF 

 function RTPDF(correctRTs, wrongRTs, missRTs, kernelRTMinMS, kernelRTMaxMS, minRespTimeMS, maxRespTimeMS)

    subplot(4, 3, 3);
    cla;
    cdfplot([correctRTs wrongRTs missRTs]);
%     timeLimit = min(file.responseLimitMS, 5000);
    timeLimit = 1000;
    set(gca, 'XLim', [-500 timeLimit], 'YLim', [0 1]);
    hold on;
    yLimits = get(gca, 'YLim');
    plot([0 0], yLimits, 'k');
    plot(double(minRespTimeMS) * [1 1], yLimits, '--', 'Color', 0.5 * [0 1 0]);
    plot(double(maxRespTimeMS) * [1 1], yLimits, '--', 'Color', 0.5 * [1 0 0]);
    plot(double(kernelRTMinMS) * [1 1], yLimits, ':', 'Color', 0.5 * [0 1 0]);
    plot(double(kernelRTMaxMS) * [1 1], yLimits, ':', 'Color', 0.5 * [1 0 0]);
    xlabel('Time Relative to Stimulus');
    ylabel('');
    title('Cumulative Reaction Times');
 end