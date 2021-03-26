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
	[U, ~] = getSubset('normal', dataDirName, tableDataName, [], limits);
%   for i = 1:height(U)
  for i = 26:100
    processFile(dataDirName + U.animal(i) + '/MatFiles/' + U.date(i) + '.mat');
  end
end

function processFile(fileName)
  load(fileName, 'file', 'trials');
  OKMatlab([], file, trials);                             % display the results
  
  indices.correct = find([trials(:).trialEnd] == 0);      % correct trials
  indices.fa = find([trials(:).trialEnd] == 1);           % false alarm trials
  indices.miss = find([trials(:).trialEnd] == 2);         % miss trials
  numFs = length(indices.fa);
  trialStructs = [trials(:).trial];                       % trial structs extracted from trials array
  preStimMS = [trialStructs(:).preStimMS];                % planned stimulus on time
  tooFastMS = file.tooFastMS;
  
  stimOnHMS = preStimMS(indices.correct);                 % stimOn for hits
  stimOnMMS = preStimMS(indices.miss);
  stimOnFMS = preStimMS(indices.fa);
  RTs = [trials(indices.fa).reactTimeMS];
  if length(RTs) > length(indices.fa)
    for t = 1:length(indices.fa)
      if length(trials(t).reactTimeMS) > 1
        trials(t).reactTimeMS = trials(t).reactTimeMS(1);
      end
    end
    RTs = [trials(indices.fa).reactTimeMS];
  end
  FTrialTimesMS = stimOnFMS + RTs;
  
%   FTrialTimeMS = FTrialTimeMS(FTrialTimeMS > file.preStimMinMS);
    
  dT = 250;
%   histEdges = 1:tInc:file.preStimMaxMS + tooFastMS;
  edges = file.preStimMinMS:dT:file.preStimMaxMS;
  % compile a histogram of the hit rate as a function of stimOn time.
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
  doOneHist(3, stimOnHMS, sprintf('Hits (n=%d)', length(indices.correct)), 'StimOn Time', hColor);
  doOneBar(4, sprintf('H+M (n=%d)', length(indices.correct) + length(indices.miss)), 'StimOn Time', file, dT, ...
    {stimOnMMS, stimOnHMS}, {mColor, hColor});
	doOneHist(5, stimOnMMS, sprintf('Misses (n=%d)', length(indices.miss)), 'StimOn Time', mColor);
  doOneBar(6, sprintf('H+M+FA (n=%d)', length(indices.correct) + length(indices.miss) + numFs), ...
    'StimOn Time (reached or unreached)', file, dT, {stimOnFMS, stimOnMMS, stimOnHMS}, {fColor, mColor, hColor});
	doOneHist(7, stimOnFMS, sprintf('FAs (n=%d)', numFs), '(Unreached) StimOn for FA', fColor);

  % Produce a linear fit to trials unavailable for false alarms because they were consumed by hits or misses. 
  % The normal assumption is that available trials go to zero at preStimMaxMS, but that is not always the case. 
  
%   binEdgesMS = file.preStimMinMS:dT:file.preStimMaxMS;
% 	x = binEdgesMS(1:end - 1) + dT / 2.0;
%   y = histcounts([stimOnHMS, stimOnMMS], binEdgesMS);   % combined hit and miss stimOn distributions
%   coeff = polyfit(x, y, 1);
%   preStimMaxFrac = ((file.preStimMaxMS + tooFastMS) * coeff(1) + coeff(2)) / coeff(2);
%   subplot(4, 2, 4);
%   hold on;
%   plot([file.preStimMinMS, file.preStimMaxMS], ...
%     [file.preStimMinMS * coeff(1) + coeff(2), file.preStimMaxMS * coeff(1) + coeff(2)]);
  
  % adjust for the abbreviated trials
%   binEdgesMS = file.preStimMinMS + tooFastMS:dT:file.preStimMaxMS + tooFastMS;
%   binEdgesMS = file.preStimMinMS:dT:file.preStimMaxMS + tooFastMS;
% 	x = binEdgesMS(1:end - 1) + dT / 2.0;
%   y = histcounts(FTrialTimesMS, binEdgesMS);
%   y = adjustForShortTrials(y, binEdgesMS, file.preStimMinMS + tooFastMS, file.preStimMaxMS + tooFastMS, preStimMaxFrac);
   [x, y, varY, n, l2] = correctedFACounts(FTrialTimesMS, [stimOnHMS, stimOnMMS], file, dT);

   % fit() doesn't work well if there are a lot of y == 0 entries at the end of the sequence.  These are undersampled,
  % and they can't be corrected for the number of short trials. 
%   lastY = find(y > 0, 1, 'last' );
%   y = y(1:lastY)';
%   x = x(1:lastY)';
  fitted = fit(x', y', 'exp1', 'startPoint', [0, 0], 'weight', n);
  subplot(4, 2, 8);
  hold on;
  plot(fitted, x, y);
  hl = findobj(gcf,'type','legend');
  delete(hl);
  xlabel('Time (ms)');
  ylabel('');
  axis([0, file.preStimMaxMS + tooFastMS, 0, inf]);
  doOneHist(8, FTrialTimesMS, sprintf('FAs (n=%d)', numFs), 'FA Time from Trial Start (ms)', fColor);
  
  hitRate = length(indices.correct) / (length(indices.correct) + length(indices.miss));
%   header(file, hitRate, -fitted.b);
  header(file, hitRate, l2 / 1000.0);

  sameYAxisScaling(4, 2, [3, 5, 7]);
  sameYAxisScaling(4, 2, [4, 6]);

end

%%
% function y = adjustForShortTrials(y, binEdgesMS, minMS, maxMS, preStimMaxFrac)
% %
% % Adjust for abbreviated trials
% %
%   % The bins that occur entirely before minMS need no adjustment at all
%   b = 1;
%   while b < length(binEdgesMS) - 1 && binEdgesMS(b + 1) < minMS 
%     b = b + 1;
%   end
%   
%   % The first bin spans minMS needs special handling.  The fraction of that bin before minMS needs no adjustment,
%   % but a proportional adjustment is needed for the fraction after minMS.
%   adjust = (minMS - binEdgesMS(b)) / diff(binEdgesMS(b:b + 1));              % portion before minMS
%   fullFrac = 1.0 - mean([binEdgesMS(b + 1), minMS])  * (1.0 - preStimMaxFrac) / (maxMS - minMS); % frac full trials in changing part of bin
%   adjust = adjust + fullFrac * (binEdgesMS(b + 1) - minMS) / diff(binEdgesMS(b:b + 1));
%   y(b) = y(b) / adjust;
%   b = b + 1;
%   
%   % The remaining bins have a linearly changing number of short trials across their full width.
%   for b = b:length(binEdgesMS) - 1
%     fullFrac = 1.0 - (mean(binEdgesMS(b:b + 1)) - minMS) * (1.0 - preStimMaxFrac) / (maxMS - minMS); % average full trials across entire bin
%     y(b) = y(b) / fullFrac;
%   end
% end

%%
function [x, y, varY, numF, l2] = correctedFACounts(FTrialTimesMS, HMStimOnTimesMS, file, dT)
%   
% The times after hit/miss truncations need to be incremented in proportion to the fraction of trials
% that have been truncated by hit/miss trials at that point. The fractional increment is based on a count
% that excludes all trial for which a FA has already occurred for abbreviated trials
%
  numTrials = length(FTrialTimesMS) + length(HMStimOnTimesMS);
	binEdgesMS = 0:dT:file.preStimMaxMS;
	x = binEdgesMS(1:end - 1) + dT / 2.0;
    
  fprintf('Total Mean rate %.2f per second\n', 1000.0 / mean(FTrialTimesMS));
  pastPreStim = FTrialTimesMS > file.preStimMinMS;
  sumTime = sum(FTrialTimesMS(pastPreStim));
  sumTime = sumTime + sum(HMStimOnTimesMS - file.preStimMinMS);
  lambdaS = sumTime / length(FTrialTimesMS) / 1000.0;
  
%   lambdaS = 1000.0 / mean(FTrialTimesMS(pastPreStim));
  fprintf('Prestim Onward Mean rate %.2f per second\n', lambdaS);
  fprintf('  Sum post prestim: %.2f\n', sumTime);
  fprintf('    N post prestim: %.0f\n', length(FTrialTimesMS));
  integral = 1.0 / lambdaS;
  fprintf('  integral  %.2f\n', integral);
  truncX = exp(-lambdaS * file.preStimMaxMS / 1000.0);
  truncIntegral = 1.0 / lambdaS * (truncX * log(truncX) - truncX + 1) + truncX * file.preStimMaxMS / 1000.0;
  l2 = lambdaS;
  while truncIntegral < integral
    l2 = l2 - 0.01;
    truncX = exp(-l2 * file.preStimMaxMS / 1000.0);
    truncIntegral = 1.0 / l2 * (truncX * log(truncX) - truncX + 1) + truncX * file.preStimMaxMS / 1000.0;
    fprintf('lambdaS %.4f, l2 %.4f, integral %.4f truncIntegral %.4f\n', lambdaS, l2, integral, truncIntegral);
  end
%   truncIntegral = 1.0 / lambdaS * (
% Create 1 ms histograms for the times of FAs and the times when stim on truncated hit/miss trials.  We add
% tooFastMS because a stimOn doesn't prevent a FA until after tooFastMS
  fHist = zeros(1, file.preStimMaxMS + file.tooFastMS);
  hmHist = zeros(1, file.preStimMaxMS + file.tooFastMS);
  bins = length(FTrialTimesMS);
  for b = 1:bins
    if FTrialTimesMS(b) <= length(fHist)
      fHist(FTrialTimesMS(b)) = fHist(FTrialTimesMS(b)) + 1;
    end
  end
  bins = length(HMStimOnTimesMS);
  for b = 1:bins
    hmHist(HMStimOnTimesMS(b) + file.tooFastMS) = hmHist(HMStimOnTimesMS(b) + file.tooFastMS) + 1;
  end
  fCumHist = cumsum(fHist);
  hmCumHist = cumsum(hmHist);
  % weighting function that corrects for all the trials eliminated by hit/miss trials
  fWeight = (numTrials - fCumHist) ./ (numTrials - fCumHist - hmCumHist);
  
  figure(3);
  clf;
  plot(numTrials - fCumHist);
  hold on;
  plot(numTrials - fCumHist - hmCumHist);
  ylim([0, inf]);
  
  figure(4);
  clf;
  plot(fWeight);
  
  corrected = fHist .* fWeight;
  
  figure(5);
  clf;
  plot(corrected);
  
  figure(2)
  
  % Bin the weighted FAs using the bin edges provided
  y = zeros(1, length(binEdgesMS) - 1);
  numF = zeros(1, length(binEdgesMS) - 1);

  for b = 1:length(binEdgesMS) - 1
    y(b) = sum(corrected(max(1, binEdgesMS(b)):min(length(corrected), binEdgesMS(b + 1) - 1)));
    numF(b) = sum(fHist(max(1, binEdgesMS(b)):min(length(corrected), binEdgesMS(b + 1) - 1)));
  end
  % fit() doesn't work well if there are a lot of y == 0 entries at the end of the sequence.  These are undersampled,
  % and they can't be corrected for the number of short trials. 
  lastY = find(y > 0, 1, 'last' );
  y = y(1:lastY);
  x = x(1:lastY);
  numF = numF(1:lastY);
  
  pFA = y / numTrials / (length(binEdgesMS) - 1);
  varY = (pFA .* (1.0 - pFA)) ./ numF;
  l2 = 1.0 / l2;
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
function [dP, c] = dPrime(pH, pFA)
%
% Compute d-prime based on hits and false alarms 
%
% convert to Z scores and calculate d-prime and criterion
  zH = -sqrt(2) .* erfcinv(2 * pH);
  zFA = -sqrt(2) .* erfcinv(2 * pFA);
  dP = zH - zFA ;
  c = -0.5 * (zH + zFA);
end

%%
function header(file, pH, lambda)

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
    headerText{length(headerText) + 1} = sprintf('Response window: %d -- %d ms', ...
        file.tooFastMS, file.rewardedLimitMS);
  elseif isfield(file, 'tooFastMS') && isfield(file, 'reactMS')
    headerText{length(headerText) + 1} = sprintf('Response window: %d -- %d ms', ...
        file.tooFastMS, file.reactMS);
  end
	headerText{length(headerText) + 1} = sprintf('FA/s %.2f', lambda * 1000.0);
	pFA = 1.0 - exp(-lambda * (file.rewardedLimitMS - file.tooFastMS));
	headerText{length(headerText) + 1} = sprintf('FA rate %.2f', pFA);
  headerText{length(headerText) + 1} = sprintf('Raw hit rate %.2f', pH);
  pHtrue = (pH - pFA) / (1.0 - pFA);
	headerText{length(headerText) + 1} = sprintf('True hit rate %.2f', pHtrue);
  [dP, c] = dPrime(pHtrue, pFA);
	headerText{length(headerText) + 1} = sprintf('d-prime = %.2f; criterion %.2f', dP, c);

  text(0.00, 1.00, headerText, 'VerticalAlignment', 'top');
end