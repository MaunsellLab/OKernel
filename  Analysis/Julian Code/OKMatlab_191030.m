function dParams = OKMatlab_191030(dParams, file, trials)

% OKMatlab is invoked at the end of every trial.  As a free-standing function that is instantiated on each call,
% it has no way to store static across calls.  Instead, such values are stored in a struct, dParams.  dParams
% arrives as an empty matrix on the first call, so the first call can be identified in this way.  By returning
% dParam with essential values, they can be recovered in each call.

h = findobj(allchild(groot), 'flat', 'type', 'figure', 'number', 1);
% Set up the figure if it hasn't been prepared yet
if isempty(h) || ~strcmp(get(h, 'type'), 'figure') || nargin == 1 
    clear dParams;
%    dParams.startTimeVec = now;
    dParams.plotLayout = {4,3};
    dParams.RTBins = 10;
    figure('Units', 'inches', 'Position', [20, 5, 8.5, 11]);
    clf;
    drawText(dParams);
    if nargin == 1
        return;
    end
end

% enforce the subjectNumber as a scalar
if ~isscalar(file.subjectNumber)
file.subjectNumber = file.subjectNumber(end);
end

file.trials = size(trials, 2);              % update the number of trials
if length([trials(:).trialEnd]) ~= length(trials)
    disp('Adding trialEnd to trial(s) without');
    for i = 1:length(trials)
        if isempty(trials(i).trialEnd)
            trials(i).trialEnd = -1;
        end
    end
end
for t = 1:length(trials)
    if length(trials(t).reactTimeMS) > 1
        trials(t).reactTimeMS = trials(t).reactTimeMS(1);
    end
end

% all the functions need to have the current trial ends

indices.correct = find([trials(:).trialEnd] == 0);      % correct trials
indices.fa = find([trials(:).trialEnd] == 1);           % false alarm trials
indices.miss = find([trials(:).trialEnd] == 2);         % miss trials
trialStructs = [trials(:).trial];                       % trial structs extracted from trials array

drawText(dParams, file, trials);
RTHistogram(dParams, file, trials, indices);
RTPDF(dParams, file, trials, indices);
cumulativeTime(dParams, trials);
outcomesOverDay(dParams, trials, indices);
outcomesOverTrial(dParams, file, trials, indices, trialStructs);

RTvStimulus(dParams, trials, trialStructs);
psychometric(dParams, trials, trialStructs);
oKernel(dParams, file, trials);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performance Text

function drawText(dParams, file, trials)

axisHandle = subplot(dParams.plotLayout{:}, 1);						% default axes are 0 to 1
set(axisHandle, 'Visible', 'off');
set(axisHandle, 'OuterPosition', [0.02 0.75, 0.25, 0.2]);
text(0.00, 1.25, 'OKernel', 'FontWeight', 'bold', 'FontSize', 16);
if (nargin == 1)                                % first call of function, display task title only
    return;
end

%elapsedMinutes = round((now - datenum(file.startTimeVec)) * 24 * 60);
%startStr = datestr(file.startTimeVec,'mmmm dd, yyyy HH:MM');
%pmStr = '\pm';
muStr = '\mul';
if file.subjectNumber == 0 || file.subjectNumber == 999
    text(0.00, 1.11, sprintf('Subject: %d', file.subjectNumber), 'FontSize', 14, 'Color', 'r');
else
    text(0.00, 1.11, sprintf('Subject: %d', file.subjectNumber), 'FontSize', 14);
end
headerText = cell(1, 1);
if isfield(file, 'startTimeVec')
    headerText{1} = datestr(file.startTimeVec, 'mmmm dd, yyyy HH:MM');
else
    headerText{1} = '(missing date field)';
end
if isfield(trials(:), 'rewardUL')
    headerText{length(headerText) + 1} = sprintf('Task active for %.0f minutes, %.0f %s', ...
        (trials(end).totalRunTimeS - trials(1).totalRunTimeS) / 60.0, sum([trials(:).rewardUL]), muStr);
else
    headerText{length(headerText) + 1} = sprintf('Task active for %.0f minutes',...
        (trials(end).totalRunTimeS - trials(1).totalRunTimeS) / 60.0);
end
if isfield(file, 'tooFastMS') && isfield(file, 'rewardedLimitMS')
    headerText{length(headerText) + 1} = sprintf('Too fast to response limit: %d -- %d ms', ...
        file.tooFastMS, file.rewardedLimitMS);
end
headerText{length(headerText) + 1} = sprintf('Pulse Width: %d ms', round(trials(1).trial.pulseDurMS));
if isfield(file, 'kernelRTMinMS') && isfield(file, 'kernelRTMaxMS')
    headerText{length(headerText) + 1} = sprintf('Kernel RT limits: %d -- %d ms', file.kernelRTMinMS, file.kernelRTMaxMS);
end
text(0.00, 1.00, headerText, 'VerticalAlignment', 'top');

correctIndices = find([trials(:).trialEnd] == 0);
numCorrect = length(correctIndices);
failedIndices = find([trials(:).trialEnd] == 2);
numFailed = length(failedIndices);
earlyIndices = find([trials(:).trialEnd] == 1);
numEarly = length(earlyIndices);
totalTrials = numCorrect + numFailed + numEarly;

t2H(1) = text(0.00, 0.55, {'Trials:', 'Correct:', 'Failed:', 'Early:'});
t2H(2) = text(0.35, 0.55, {sprintf('%d', file.trials), sprintf('%d', numCorrect), sprintf('%d', numFailed), ...
                    sprintf('%d', numEarly)});
t2H(3) = text(0.60, 0.55, {' ', sprintf('%.0f%%', numCorrect / totalTrials * 100.0), ...
        sprintf('%.0f%%', numFailed / totalTrials * 100.0), sprintf('%.0f%%', numEarly / totalTrials * 100.0)});
set(t2H, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
set(gcf, 'Visible', 'on');

if trials(end).trial.syntheticData == 1
    text(0.00, 0.10, 'Synthetic Data', 'FontSize', 14, 'Color', 'r', 'FontWeight', 'bold', 'VerticalAlignment', 'top');
end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% RT Histogram %%%

 function RTHistogram(dParams, file, trials, indices)
 if isempty(indices.correct)
     correctRTs = -10000;                  % make sure we don't get an empty matrix from histc
 else
     correctRTs = [trials(indices.correct).reactTimeMS];
 end
 if isempty(indices.fa)
    wrongRTs = -10000;                  % make sure we don't get an empty matrix from histc
 else
    wrongRTs = [trials(indices.fa).reactTimeMS];
 end
 if isempty(indices.miss)
    missRTs = -10000;                  % make sure we don't get an empty matrix from histc
 else
    allMissRTs = [trials(indices.miss).reactTimeMS];
    missRTs = allMissRTs(allMissRTs > 0);
 end
 subplot(dParams.plotLayout{:}, 7);
 cla;
 hold on;
 timeLimit = min(file.responseLimitMS, 5000);
 edges = linspace(-1000, timeLimit, dParams.RTBins);
 nCorrect = histc(correctRTs, edges);
 nWrong = histc(wrongRTs, edges);
 nMiss = histc(missRTs, edges);
 if sum(nCorrect) + sum(nWrong) + sum(nMiss) > 0
    binSize = edges(2) - edges(1);
    bH = bar(edges + binSize / 2, [nCorrect(:), nMiss(:), nWrong(:)], 'stacked');
    set(bH, 'barWidth', 1, 'lineStyle', 'none');
    set(bH(1), 'FaceColor', [0 0 0.6]);
    set(bH(2), 'FaceColor', [0.6 0 0]);
    set(bH(3), 'faceColor', [0.6 0 0]);
    if max([nWrong, nCorrect, nMiss] > 50)         % re-bin on next plot?
       dParams.RTBins = min([dParams.RTBins * 2, 100]);
    end
    yLimits = get(gca, 'YLim');                % vertical line at stimulus on
    plot([0 0], yLimits, 'k');
    if isfield(file, 'tooFastMS')
        llH = plot(double(file.tooFastMS) * [1 1], yLimits, 'k--');
        set(llH, 'Color', 0.5 * [0 1 0]);
    end
    if isfield(file, 'rewardedLimitMS')
        llH = plot(double(file.rewardedLimitMS) * [1 1], yLimits, 'r--');
        set(llH, 'Color', 0.5 * [1 0 0]);
    end
    if isfield(file, 'responseLimitMS')
        llH = plot(double(file.responseLimitMS) * [1 1], yLimits, 'r--');
        set(llH, 'Color', 0.5 * [1 0 0]);
    end
 end
 set(gca, 'XLim', [-1000 timeLimit]);
 xlabel('Time Relative to Stimulus');
 title('Reaction Times');

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% RT PDF %%%

 function RTPDF(dParams, file, trials, indices)

 if isempty(indices.correct)
     correctRTs = -10000;                  % make sure we don't get an empty matrix from histc'
 else
     correctRTs = [trials(indices.correct).reactTimeMS];
 end
 if isempty(indices.fa)
    wrongRTs = -10000;                    % make sure we don't get an empty matrix from histc
 else
    wrongRTs = [trials(indices.fa).reactTimeMS];
 end
 if isempty(indices.miss)
    missRTs = -10000;                    % make sure we don't get an empty matrix from histc
 else
    missRTs = [trials(indices.miss).reactTimeMS];
%     missRT(missRTs < 0) = 100000;        % include misses in count, but don't let them display on plot
 end
 subplot(dParams.plotLayout{:}, 4);
 cla;
 cdfplot([correctRTs wrongRTs missRTs]);
 timeLimit = min(file.responseLimitMS, 5000);
 set(gca, 'XLim', [-1000 timeLimit], 'YLim', [0 1]);
 hold on;
 yLimits = get(gca, 'YLim');
 plot([0 0], yLimits, 'k');
 if isfield(file, 'tooFastMS')
    plot(double(file.tooFastMS) * [1 1], yLimits, '--', 'Color', 0.5 * [0 1 0]);
 end
 if isfield(file, 'rewardedLimitMS')
    plot(double(file.rewardedLimitMS) * [1 1], yLimits, '--', 'Color', 0.5 * [1 0 0]);
 end
 if isfield(file, 'responseLimitMS')
    plot(double(file.responseLimitMS) * [1 1], yLimits, '--', 'Color', 0.5 * [1 0 0]);
 end
 if isfield(file, 'kernelRTMinMS')
    plot(double(file.kernelRTMinMS) * [1 1], yLimits, ':', 'Color', 0.5 * [0 1 0]);
 end
 if isfield(file, 'kernelRTMaxMS')
    plot(double(file.kernelRTMaxMS) * [1 1], yLimits, ':', 'Color', 0.5 * [1 0 0]);
 end

 title('Cumulative Reaction Times');
 xlabel('');
 ylabel('');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Display the optogenetic kernel %%%

function oKernel(~, file, trials)
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

function hasSubplot = figureHasSubplotTitleString(figureNum, title)
    
    h = findobj(allchild(groot), 'flat', 'type', 'figure', 'number', figureNum);
    hasSubplot = false;
    for s = 1:length(h.Children)
        if ~isempty(strfind(h.Children(s).Title.String, title))
            hasSubplot = true;
            break
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RT versus Stimulus Intensity %%%

function RTvStimulus(dParams, trials, trialStructs)

hit = [trials(:).trialEnd] == 0;
stimValueSet = unique([trials(hit).visualStimValue]);
numStim = length(stimValueSet);
noStim =  [trialStructs(:).optoIndex] == 0;
delay0 = [trialStructs(:).optoIndex] == 1;
delay1 = [trialStructs(:).optoIndex] == 2;
noStimRTMean = zeros(1, numStim);
noStimRTSE = zeros(1, numStim);
delay0RTMean = zeros(1, numStim);
delay0RTSE = zeros(1, numStim);
delay1RTMean = zeros(1, numStim);
delay1RTSE = zeros(1, numStim);
for s = 1:length(stimValueSet)
    stimTrial = [trials(:).visualStimValue] == stimValueSet(s);
    noStimRTMean(s) = mean([trials(stimTrial & hit & noStim).reactTimeMS]);
    noStimRTSE(s) = std([trials(stimTrial & hit & noStim).reactTimeMS]) / sqrt(sum(stimTrial & hit & noStim));
    delay0RTMean(s) = mean([trials(stimTrial & hit & delay0).reactTimeMS]);
    delay0RTSE(s) = std([trials(stimTrial & hit & delay0).reactTimeMS]) / sqrt(sum(stimTrial & hit & delay0));
    delay1RTMean(s) = mean([trials(stimTrial & hit & delay1).reactTimeMS]);
    delay1RTSE(s) = std([trials(stimTrial & hit & delay1).reactTimeMS]) / sqrt(sum(stimTrial & hit & delay1));
end
subplot(dParams.plotLayout{:}, 8);
cla;
hold off;
errorbar(stimValueSet, delay0RTMean, delay0RTSE, '-s', 'color', [0.8500, 0.3250, 0.0980], 'markersize', 6, ...
         'markerfacecolor', [0.8500, 0.3250, 0.0980], 'markeredgecolor', [0.8500, 0.3250, 0.0980]);
set(gca,'xscale', 'log', 'xgrid', 'on');
hold on;
errorbar(stimValueSet, delay1RTMean, delay1RTSE, '-s', 'color', [0.6350, 0.0780, 0.1840], 'markersize', 6, ...
         'markerfacecolor', [0.6350, 0.0780, 0.1840], 'markeredgecolor', [0.6350, 0.0780, 0.1840]);
hold on;
errorbar(stimValueSet, noStimRTMean, noStimRTSE, '-sb', 'markersize', 6, 'markerfacecolor', 'blue', ...
          'markeredgecolor', 'blue');
set(gca, 'XTickLabel', num2str(transpose(get(gca, 'XTick'))))          % get rid of scientific notation
yLimits = get(gca, 'YLim');
set(gca, 'YLim', [0 yLimits(2) * 1.05]);
title('Mean Reaction Times');
ylabel('RT');
if trials(end).blockStatus.visualStimType == 0
    xlabel('Contrast (%)');
else
    xlabel('Power (mW)');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Psychometric Function %%%

function psychometric(dParams, trials, trialStructs)

hit = [trials(:).trialEnd] == 0;
miss = [trials(:).trialEnd] == 2;
stimValueSet = unique([trials(hit | miss).visualStimValue]);
numStim = length(stimValueSet);
delay0 = [trialStructs(:).optoIndex] == 1;
delay0Hits = zeros(1, numStim);
delay0N = zeros(1, numStim);
%delay1 = [trialStructs(:).optoIndex] == 2;
%delay1Hits = zeros(1, numStim);
%delay1N = zeros(1, numStim);
noStimHits = zeros(1, numStim);
noStimN = zeros(1, numStim);
for s = 1:numStim                                          % for each stim value
    stimTrial = [trials(:).visualStimValue] == stimValueSet(s);
    delay0Hits(s) = sum(stimTrial & hit & delay0);
    delay0N(s) = delay0Hits(s) + sum(stimTrial & miss & delay0);
    noStimHits(s) = sum(stimTrial & hit & ~delay0);
    noStimN(s) = noStimHits(s) + sum(stimTrial & miss & ~delay0);
end
if ~isempty(delay0Hits)
    [delay0HitRate, delay0pci] = binofit(delay0Hits, delay0N);
    delay0YNeg = delay0HitRate - delay0pci(:, 1)';
    delay0YPos = delay0pci(:, 2)' - delay0HitRate;
end
if ~isempty(noStimHits)
    [noStimHitRate, noStimPci] = binofit(noStimHits, noStimN);
    noStimYNeg = noStimHitRate - noStimPci(:, 1)';
    noStimYPos = noStimPci(:, 2)' - noStimHitRate;
end
%delay1YNeg = delay1HitRate - delay1pci(:, 1)';
%delay1YPos = delay1pci(:, 2)' - delay1HitRate;
misses = length(miss);
missRate = misses / (misses + length(hit));
falseHitRate = trials(end).falseHitRate / (trials(end).falseHitRate + missRate);

subplot(dParams.plotLayout{:}, 5);
cla;
if exist('noStimHitRate', 'var') == 1
    hold off;
    errorbar(stimValueSet, noStimHitRate, noStimYNeg, noStimYPos, '-s', 'markersize', 6, 'markerfacecolor', 'blue');
    hold on;
end
if exist('delay0HitRate', 'var') == 1
    errorbar(stimValueSet, delay0HitRate, delay0YNeg, delay0YPos, '-s', 'markersize', 6, ...
         'color', [0.8500, 0.3250, 0.0980], 'markerfacecolor', [0.8500, 0.3250, 0.0980]);
    hold on;
end
set(gca,'xscale','log', 'xgrid', 'on');
hold on;
set(gca, 'ylim', [0 1]);
xLimits = get(gca, 'XLim');
plot([xLimits(1), xLimits(2)], [falseHitRate, falseHitRate], 'k:');
set(gca, 'XLim', [xLimits(1), min(xLimits(2), 100)]);
set(gca, 'XTickLabel',num2str(transpose(get(gca, 'XTick'))))          % get rid of scientific notation
ylabel('Percent Correct');
title(sprintf('Hit Rates, 95%% CI (n = %d - %d)', min([delay0N, noStimN]), max([delay0N, noStimN])));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Work Time Efficiency %%%

function cumulativeTime(dParams, trials)

minRunTimeS = [trials(:).minRunTimeS] - trials(1).minRunTimeS(1);
totalRunTimeS = [trials(:).totalRunTimeS] - trials(1).totalRunTimeS(1);
axisHandle = subplot(dParams.plotLayout{:}, 2);
cla;
hold off;
plot(totalRunTimeS / 60.0, minRunTimeS / 60.0, 'b-', 'linewidth', 1.0);
hold on;
xlabel('Actual work time (min)');
ylabel('Ideal work time (min)');
a = axis;
limit = max(a(2), a(4));
set(axisHandle, 'xlim', [0, limit], 'ylim', [0, limit]);
plot([0, limit], [0, limit], 'k--');
efficiency =  mean(minRunTimeS) / mean(totalRunTimeS);
plot([0, limit], [0, limit * efficiency], 'k:');
title(sprintf('Work Efficiency %.0f%%', efficiency * 100.0));

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% Trial Outcomes Within Trials %%%

 function outcomesOverTrial(dParams, file, trials, indices, trialStructs)

 preStimMS = [trialStructs(:).preStimMS];
 hitTimes = preStimMS(indices.correct);
 hitRTs = [trials(indices.correct).reactTimeMS];    % release relative to stimTime on hits
 faTimes = preStimMS(indices.fa);
 faRTs = [trials(indices.fa).reactTimeMS];          % response relative to stimTime on FAs
 missTimes = preStimMS(indices.miss);
 if isempty(hitTimes)
    hitTimes = -10000;
    hitRTs = 0;
 end
 if isempty(faTimes)
     faTimes = -10000;
     faRTs = 0;
 end
 if isempty(missTimes)
     missTimes = -10000;
 end
 releaseTimes = [(hitTimes + hitRTs), (faTimes + faRTs)];

 % histc uses edges in a ridiculous way.  The first and last edge limit the range, but the bin edges used in between
 % aren't the remaining values.  They are the spots halfway between the remaining values.  Absurd.  This is presumably
 % why Matlab recommends histcounts.  But we don't have that function on all our versions of Matlab.  We set the bin
 % edges up so that there are 2 small bins on either side of the range, which we will merge at the end of the counting.

 trialBins = max(10, length(trials) / 20);                    	% 10 or more bins
 timeRangeMS = file.preStimMaxMS + file.responseLimitMS;
 binWidth = timeRangeMS / trialBins;                            % binWidth in MS
 edges = [0.0, linspace(0.5 * binWidth, timeRangeMS - 0.5 * binWidth, trialBins), timeRangeMS];
 counts = hist(hitTimes, edges);
 nHits = [counts(1) + counts(2), counts(3:end-2), counts(end-1) + counts(end)];
 counts = hist(faTimes, edges);
 nFAs = [counts(1) + counts(2), counts(3:end-2), counts(end-1) + counts(end)];
 counts = hist(missTimes, edges);
 nMisses = [counts(1) + counts(2), counts(3:end-2), counts(end-1) + counts(end)];
 nTotal = nHits + nFAs + nMisses;
 counts = hist(releaseTimes, edges);
 nReleases = [counts(1) + counts(2), counts(3:end-2), counts(end-1) + counts(end)];
 if sum(nTotal) < 1                                     % no counts yet?
     return;
 end

 pHits = nHits ./ nTotal;                               % proportions of hits, FAs and misses
 pFAs = nFAs ./ nTotal;
 pMisses = nMisses ./ nTotal;
 pReleases = nReleases ./ (max(nReleases) * 1.25);
 subplot(dParams.plotLayout{:}, 6);
 cla;
 hold on;
 set(gca, 'XLim', [1 trialBins]);
 minStimXPos = file.preStimMinMS / timeRangeMS * trialBins + 1;
 maxStimXPos = file.preStimMaxMS / timeRangeMS * trialBins + 1;
 set(gca,'XTick', [1, minStimXPos, maxStimXPos, trialBins]);
 set(gca,'XTickLabel',{'0', sprintf('%d', file.preStimMinMS), sprintf('%d', file.preStimMaxMS), ...
     sprintf('%d', timeRangeMS)});
 xlabel('Time of stimulus onset (ms)');
 ylabel('Proportion of trials');
 set(gca, 'YLim', [0 1]);
 set(gca,'YTick', [0, 1]);
 title('Outcomes Over Trial');

 plot(pHits, 'color', [0.0, 0.7, 0.0], 'lineWidth', 1);
 plot(pFAs, 'color', [0.9, 0.0, 0.0], 'lineWidth', 1);
 plot(pMisses, 'color', [0.6, 0.4, 0.2], 'lineWidth', 1);
 plot(pReleases, 'color', [0.6, 0.6, 0.6], 'lineWidth', 1);
 plot([minStimXPos, minStimXPos], [0, 1], 'k:');
 plot([maxStimXPos, maxStimXPos], [0, 1], 'k:');
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% Trial Outcomes Over Day %%%
 
 function outcomesOverDay(dParams, trials, indices)
 
 numTrials = length(trials);
 hits = zeros(1, numTrials);
 fas = zeros(1, numTrials);
 misses = zeros(1, numTrials);
 hits(indices.correct) = 1;
 fas(indices.fa) = 1;
 misses(indices.miss) = 1;
 
 subplot(dParams.plotLayout{:}, 3);
 cla;
 hold on;
 plot(smooth(hits, ceil(numTrials / 5), 'lowess'), 'color', [0.0, 0.7, 0.0], 'lineWidth', 1);
 plot(smooth(fas, ceil(numTrials / 5), 'lowess'), 'color', [0.9, 0.0, 0.0], 'lineWidth', 1);
 plot(smooth(misses, ceil(numTrials / 5), 'lowess'), 'color', [0.6, 0.4, 0.2], 'lineWidth', 1);
 
 xlabel('Trials');
 ylabel('Proportion of trials');
 set(gca, 'YLim', [0 1]);
 set(gca,'YTick', [0, 1]);
 title('Outcomes Over Day');
