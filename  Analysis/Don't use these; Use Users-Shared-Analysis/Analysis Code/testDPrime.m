function testDPrime()

%{
It's not straightforward to get a d-Prime from this task design.  We can't simply related the H to
the number of FAs.  An animals with a very large d' could have many more FAs than Hs if the preStim period
is very long.  Conversely, and animal with a very small d' could have many more Hs than FAs if the preStim
period is very short.  

We need to compute d' in way that takes the preStim interval into account.  This can be done by dividing
the trials into preStim bins.  It's straightforward to the number of Hs and Ms for each bin.  The remaining
(executed) trials are FAs, but they include early FAs that are not relevant to determining the probability of 
a H or FA once the animal has gotten to a give preStim bin, which is what we want. We need to take the FA
trials that reach the bin in question, and divide them between FAs that occur during the bin, and other FAs, 
which are the CR for the bin.  

We don't need to concern ourselves with the reactTimeMS.  The question we are addressing is whether the animal
got a H, M or FA on the trials when it got to the bin that contained the stimulusOn. 

A final complication is the tooFastTime.  If tooFastTime is zero, then the condition for the FA is reaching the
preStim time.  But we want to take the tooFastTime into account, so we consider only FAs that reach preStim time 
plus tooFastTime and ask what fraction of them ended in the FA before the reaction time had expired.  These 
define the FA rate that we want to compare with the H rate for the bin.

We collect the FA rates based on the preStim for each H/M trial. The resulting values (H/M/FA) are placed into fixed 
preStim bins.

%}

% Get a list of files to examine
  dataDirName = '/Users/Shared/Data/OKernel/';
	tableDataName = [dataDirName ' Analysis/Processed Files.mat'];
  limits.rampMS = 0;
  limits.criterion = 0;
  limits.animal = {'All'};
  limits.minTrials = 0;
  limits.minDec = -1;
  limits.oneDay = [];
  limits.minSessions = 0;
% 	[U, ~] = getSubset('normal', dataDirName, tableDataName, limits);
	[U, dataDirName] = getSessionTable('example');

  for i = 1:height(U)
    matFileName = dataDirName + U.animal(i) + '/MatFiles/' + U.date(i) + '.mat';
    processFile(i, U, matFileName);
    fprintf('File %d of %d pHit %.2f, pFA %.2f, d'' %.2f\n', i, height(U), U.pHit(i), U.pFA(i), U.dPrime(i));
    folderName = strcat(dataDirName, ' Analysis/Figures/D-Primes/') + U.animal(i);
    if ~exist(folderName, 'dir')
       mkdir(folderName);
    end
    saveas(gcf, folderName + '/' + U.date(i) + '.pdf')
  end
end

function processFile(i, U, matFileName)
  load(matFileName, 'file', 'trials');
%   OKMatlab([], file, trials);                             % display the results
  
  indices.correct = [trials(:).trialEnd] == 0;              % correct trials
  indices.early = [trials(:).trialEnd] == 1;                   % false alarm trials
  indices.fail = [trials(:).trialEnd] == 2;                 % miss trials
  [respLimitsMS, indices, fitCum, endCumTimeMS] = getResponseLimits(file, trials, indices);
  numFs = length(indices.early);
  trialStructs = [trials(:).trial];                         % trial structs extracted from trials array
  preStimMS = [trialStructs(:).preStimMS];                  % planned stimulus on time
  
  stimOnHMS = preStimMS(indices.correct);                   % stimOn for hits
  stimOnMMS = preStimMS(indices.fail);
  stimOnFMS = preStimMS(indices.early);
  RTs = oneRTPerTrial(trials, indices.early);
  
  FTrialTimesMS = stimOnFMS + RTs;
  dT = 250;
  edges = file.preStimMinMS:dT:file.preStimMaxMS;
  h = zeros(1, length(edges));
  n = zeros(1, length(edges));
  for t = 1:length(edges)
    h(t) = sum(stimOnHMS >= edges(t) & stimOnHMS < edges(t) + dT);
    n(t) = h(t) + sum(stimOnMMS >= edges(t) & stimOnMMS < edges(t) + dT);
  end
  if ~isempty(h)
      [pH, hPCI] = binofit(h, n);
      negCI = pH - hPCI(:, 1)';
      posCI = hPCI(:, 2)' - pH;
  end
  
  figure(2);
	set(gcf, 'units', 'inches', 'position', [27, 10.0, 7.5, 10]);    
  clf;
  hColor = [0.0, 0.7, 0.0];
  mColor = [0.3, 0.2, 0.0];
  fColor = [0.8, 0.0, 0.0];
 
  doOnePlot(2, edges, pH, negCI, posCI, 'Hit Rate', 'StimOn Time', hColor);  
	RTPDF(3, file, trials, indices, respLimitsMS, fitCum, endCumTimeMS);
  doOneBar(4, sprintf('H+M (n=%d)', sum(indices.correct) + sum(indices.fail)), 'StimOn Time', file, dT, ...
    {stimOnMMS, stimOnHMS}, {mColor, hColor});
  RTHistogram(5, file, trials, indices, endCumTimeMS);
  doOneBar(6, sprintf('H+M+FA (n=%d)', sum(indices.correct) + sum(indices.fail) + numFs), ...
    'StimOn Time (reached or unreached)', file, dT, {stimOnFMS, stimOnMMS, stimOnHMS}, {fColor, mColor, hColor});
  doOneHist(8, FTrialTimesMS, sprintf('FAs (n=%d)', numFs), 'FA Time from Trial Start (ms)', fColor);
 
  header(i, U, file);
  sameYAxisScaling(4, 2, [4, 6]);
end

%%
function doOneBar(plotIndex, theTitle, theLabel, file, dT, data, colors)
  
  subplot(4, 2, plotIndex);
 	binEdgesMS = file.preStimMinMS:dT:file.preStimMaxMS;
  barData = zeros(length(data), length(binEdgesMS) - 1);
  for d = 1:length(data)
    barData(d,:) = histcounts(data{d}, binEdgesMS)';
  end
  b = bar(binEdgesMS(1:end - 1) + dT / 2.0, barData, 'stacked', 'faceColor', 'flat', 'barWidth', 1);
  for d = 1:length(data)
    b(d).CData(:, 1:3) = repmat(colors{d}, length(binEdgesMS) - 1, 1);
  end
  a = b(1).Parent;
  a.XLim = [0, file.preStimMaxMS];
  a.XTick = 0:1000:file.preStimMaxMS;
  title(theTitle);
  xlabel(theLabel);
  xlim([0, inf]);
end

%%
function doOneHist(plotIndex, values, theTitle, theXLabel, faceColor)

  subplot(4, 2, plotIndex);
  h = histogram(values, 'binwidth', 250);
  h.FaceColor = faceColor;
  title(theTitle);
  xlabel(theXLabel);
  xlim([0, inf]);
end

function doOnePlot(plotIndex, x, y, negCI, posCI, theTitle, theXLabel, faceColor)

  subplot(4, 2, plotIndex);
	errorbar(x, y, negCI, posCI, '-s', 'markersize', 6, 'color', [0, 0, 0], 'markerfacecolor', faceColor);
  title(theTitle);
  xlabel(theXLabel);
  xlim([0, x(end)]);
  ylim([0, 1.0]);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RT PDF %%%

function RTPDF(plotIndex, file, trials, indices, respLimitsMS, fitCum, endCumTimeMS)
 
  allRTs = [[trials(indices.correct).reactTimeMS], [trials(indices.early).reactTimeMS],  [trials(indices.fail).reactTimeMS]];
  if isempty(allRTs)
    return
  end
  % Set the time limits
  startTime = -file.preStimMinMS;
  endTime = min(file.responseLimitMS, endCumTimeMS);
  
  subplot(4, 2, plotIndex);
  cla;
  cdfplot(allRTs);
  set(gca, 'XLim', [startTime, endTime], 'YLim', [0 1]);
  yLimits = get(gca, 'YLim');
  hold on;
  plot(startTime:endTime - 1, fitCum, 'Color', 0.5 * [1 0 0]);
  plot(respLimitsMS(1) * [1 1], yLimits, ':', 'Color', 0.5 * [1 0 0]);
  plot(respLimitsMS(2) * [1 1], yLimits, ':', 'Color', 0.5 * [1 0 0]);
   plot([0 0], yLimits, 'k');
  if isfield(file, 'tooFastMS')
    plot(double(file.tooFastMS) * [1 1], yLimits, '--', 'Color', 0.5 * [0 1 0]);
  end
 if isfield(file, 'rewardedLimitMS')
    plot(double(file.rewardedLimitMS) * [1 1], yLimits, '--', 'Color', 0.5 * [1 0 0]);
  elseif isfield(file, 'reactMS')
    plot(double(file.reactMS) * [1 1], yLimits, '--', 'Color', 0.5 * [1 0 0]);
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
  text(0.05, 0.95, sprintf('Resp. limits\n%.0f -- %.0f ms', respLimitsMS), 'units', 'normalized', ...
    'horizontalAlignment', 'left', 'verticalAlignment', 'top');
  title('Cumulative Reaction Times');
  xlabel('');
  ylabel('');
end

%%
function header(i, U, file)

  axisHandle = subplot(4, 2, 1);
  set(axisHandle, 'Visible', 'off');
  set(axisHandle, 'OuterPosition', [0.02 0.75, 0.25, 0.2]);
  text(0.00, 1.25, 'OKernel', 'FontWeight', 'bold', 'FontSize', 16);
  text(0.00, 1.11, sprintf('Subject: %d', file.subjectNumber), 'FontSize', 14);
  headerText = cell(1, 1);
  if isfield(file, 'startTimeVec')
      headerText{1} = datestr(file.startTimeVec, 'mmmm dd, yyyy HH:MM');
  else
      headerText{1} = '(missing date field)';
  end
  if isfield(file, 'tooFastMS') && isfield(file, 'rewardedLimitMS')
    headerText{length(headerText) + 1} = sprintf('React window: %d -- %d ms', ...
        file.tooFastMS, file.rewardedLimitMS);
  elseif isfield(file, 'tooFastMS') && isfield(file, 'reactMS')
    headerText{length(headerText) + 1} = sprintf('React window: %d -- %d ms', ...
        file.tooFastMS, file.reactMS);
  end
  headerText{length(headerText) + 1} = sprintf('Response window %.0f ms', U.RTWindowMS(i));
  headerText{length(headerText) + 1} = sprintf('pFA %.2f', U.pFA(i));
  headerText{length(headerText) + 1} = sprintf('Nostim: pHit %.2f; d'' %.2f, c %.2f', ...
    U.noStimPHit(i), U.noStimDPrime(i), U.noStimC(i));
  headerText{length(headerText) + 1} = sprintf('Stim:   pHit %.2f; d'' %.2f, c %.2f', ...
    U.stimPHit(i), U.stimDPrime(i), U.stimC(i));
	headerText{length(headerText) + 1} = sprintf('  Stim d'' change: %.2f', ...
    U.noStimDPrime(i) - U.stimDPrime(i));
  text(0.00, 1.00, headerText, 'VerticalAlignment', 'top');
end

%%
function RTHistogram(plotIndex, file, trials, indices, endCumTimeMS)

 subplot(4, 2, plotIndex);
 cla;
 hold on;
 if sum(indices.correct) == 0
     correctRTs = -10000;                  % make sure we don't get an empty matrix from histc
 else
     correctRTs = [trials(indices.correct).reactTimeMS];
 end
 if sum(indices.early) == 0
    wrongRTs = -10000;                  % make sure we don't get an empty matrix from histc
 else
    wrongRTs = [trials(indices.early).reactTimeMS];
 end
 if sum(indices.fail) == 0
    missRTs = -10000;                  % make sure we don't get an empty matrix from histc
 else
    allMissRTs = [trials(indices.fail).reactTimeMS];
    missRTs = allMissRTs(allMissRTs > 0);
 end
 timeLimit = min(file.responseLimitMS, endCumTimeMS);
 edges = linspace(-1000, timeLimit, 75);
 nCorrect = histc(correctRTs, edges); %#ok<*HISTC>
 nWrong = histc(wrongRTs, edges);
 nMiss = histc(missRTs, edges);
 if sum(nCorrect) + sum(nWrong) + sum(nMiss) > 0
    binSize = edges(2) - edges(1);
    bH = bar(edges + binSize / 2, [nCorrect(:), nMiss(:), nWrong(:)], 'stacked');
    set(bH, 'barWidth', 1, 'lineStyle', 'none');
    set(bH(1), 'FaceColor', [0 0 0.6]);
    set(bH(2), 'FaceColor', [0.6 0 0]);
    set(bH(3), 'faceColor', [0.6 0 0]);
    yLimits = get(gca, 'YLim');                % vertical line at stimulus on
    plot([0 0], yLimits, 'k');
    if isfield(file, 'tooFastMS')
        llH = plot(double(file.tooFastMS) * [1 1], yLimits, 'k--');
        set(llH, 'Color', 0.5 * [0 1 0]);
    end
    if isfield(file, 'rewardedLimitMS')
      llH = plot(double(file.rewardedLimitMS) * [1 1], yLimits, 'r--');
      set(llH, 'Color', 0.5 * [1 0 0]);
    elseif isfield(file, 'reactMS')
      llH = plot(double(file.reactMS) * [1 1], yLimits, 'r--');
      set(llH, 'Color', 0.5 * [1 0 0]);
    end
    if isfield(file, 'responseLimitMS')
        llH = plot(double(file.responseLimitMS) * [1 1], yLimits, 'r--');
        set(llH, 'Color', 0.5 * [1 0 0]);
    end
 end
 set(gca, 'XLim', [-file.preStimMinMS, timeLimit]);
 xlabel('Time Relative to Stimulus');
 title('Reaction Times');
end