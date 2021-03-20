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

  load('/Users/Shared/Data/OKernel/902/MatFiles/2019-09-26.mat', 'file', 'trials');
  OKMatlab([], file, trials);                             % display the results
  
  indices.correct = find([trials(:).trialEnd] == 0);      % correct trials
  indices.fa = find([trials(:).trialEnd] == 1);           % false alarm trials
  indices.miss = find([trials(:).trialEnd] == 2);         % miss trials
  numFAs = length(indices.fa);
  trialStructs = [trials(:).trial];                       % trial structs extracted from trials array
  preStimMS = [trialStructs(:).preStimMS];                % planned stimulus on time
%   tooFastMS = file.tooFastMS;
  
  stimOnHMS = preStimMS(indices.correct);                 % stimOn for hits
  stimOnMMS = preStimMS(indices.miss);
  stimOnFAMS = preStimMS(indices.fa);
  FATrialTimeMS = stimOnFAMS + [trials(indices.fa).reactTimeMS];
    
  tInc = 200;
  histEdges = 1:tInc:file.preStimMaxMS;
  edges = file.preStimMinMS:tInc:file.preStimMaxMS;
  
  n = histcounts(FATrialTimeMS, [histEdges file.preStimMaxMS]);
	FANumTrials = zeros(1, length(histEdges));
  for t = 1:length(edges)
    FANumTrials(t) = sum(n(1:end - t + 1));
  end
  
  % compile a histogram of the hit rate as a function of stimOn time.
  pHit = zeros(1, length(edges));
  for t = 1:length(edges)
    h = sum(stimOnHMS >= edges(t) & stimOnHMS < edges(t) + tInc);
    m = sum(stimOnMMS >= edges(t) & stimOnMMS < edges(t) + tInc);
    pHit(t) = h / (h + m);
  end
  
  figure(2);
	set(gcf, 'units', 'inches', 'position', [27, 10.0, 7.5, 10]);    
  clf;
  
  HColor = [0.0, 0.7, 0.0];
  MColor = [0.3, 0.2, 0.0];
  FColor = [0.8, 0.0, 0.0];
	dT = 250;
 
  doOnePlot(2, edges, pHit, 'o', 'Hits/(Hits+Misses)', 'Time After StimOn');  
  doOneHist(3, stimOnHMS, sprintf('Hits (n=%d)', length(indices.correct)), 'StimOn Time', HColor);
  doOneBar(4, sprintf('H+M (n=%d)', length(indices.correct) + length(indices.miss)), 'StimOn Time', file, dT, ...
    {stimOnMMS, stimOnHMS}, {MColor, HColor});
	doOneHist(5, stimOnMMS, sprintf('Misses (n=%d)', length(indices.miss)), 'StimOn Time', MColor);
  doOneBar(6, sprintf('H+M+FA (n=%d)', length(indices.correct) + length(indices.miss) + numFAs), ...
    'StimOn Time (reached or unreached)', file, dT, {stimOnFAMS, stimOnMMS, stimOnHMS}, {FColor, MColor, HColor});
	doOneHist(7, stimOnFAMS, sprintf('FAs (n=%d)', numFAs), '(Unreached) StimOn for FA', FColor);

  % Produce a linear fit to trials unavailable for false alarms because they were consumed by hits or misses. 
  % The normal assumption is that available trials go to zero at preStimMaxMS, but that is not always the case.  
  binEdgesMS = file.preStimMinMS:dT:file.preStimMaxMS;
	x = binEdgesMS(1:end - 1) + dT / 2.0;
  y = histcounts([stimOnHMS, stimOnMMS], binEdgesMS);   % combined hit and miss stimOn distributions
  coeff = polyfit(x, y, 1);
  preStimMaxFrac = (file.preStimMaxMS * coeff(1) + coeff(2)) / coeff(2);
  subplot(4, 2, 4);
  hold on;
  plot([file.preStimMinMS, file.preStimMaxMS], [(file.preStimMinMS * coeff(1) + coeff(2)), (file.preStimMaxMS * coeff(1) + coeff(2))]);
  
  % adjust for the abbreviated trials
  binEdgesMS = 0:dT:file.preStimMaxMS;
	x = binEdgesMS(1:end - 1) + dT / 2.0;
  y = histcounts(FATrialTimeMS, binEdgesMS);
  y = adjustForShortTrials(y, binEdgesMS, file.preStimMinMS, file.preStimMaxMS, preStimMaxFrac);
  
  % fit() doesn't work well if there are a lot of y == 0 entries at the end of the sequence.  These are undersampled,
  % and they can't be corrected for the number of short trials. 
  lastY = find(y > 0, 1, 'last' );
  y = y(1:lastY)';
  x = x(1:lastY)';
  fitted = fit(x, y, 'exp1', 'startPoint', [0, 0]);
  subplot(4, 2, 8);
  hold on;
  plot(fitted, x, y);
  hl = findobj(gcf,'type','legend');
  delete(hl);
  xlabel('Time (ms)');
  ylabel('');
  axis([0, file.preStimMaxMS, 0, inf]);
  text(0.95, 0.95, sprintf('L = %.5f', -fitted.b), 'units', 'normalized', 'horizontalAlignment', 'right',...
    'verticalAlignment', 'top', 'fontSize', 12);
  doOneHist(8, FATrialTimeMS, sprintf('FAs (n=%d)', numFAs), 'FA Time from Trial Start (ms)', FColor);
  
  hitRate = length(indices.correct) / (length(indices.correct) + length(indices.miss));
  header(file, hitRate, -fitted.b);
end

%%
%%
function y = adjustForShortTrials(y, binEdgesMS, minMS, maxMS, preStimMaxFrac)
%
% Adjust for abbreviated trials
%
  % The bins that occur entirely before minMS need no adjustment at all
  b = 1;
  while b < length(binEdgesMS) - 1 && binEdgesMS(b + 1) < minMS 
    b = b + 1;
  end
  
  % The first bin spans minMS needs special handling.  The fraction of that bin before minMS needs no adjustment,
  % but a proportional adjustment is needed for the fraction after minMS.
  adjust = (minMS - binEdgesMS(b)) / diff(binEdgesMS(b:b + 1));              % portion before minMS
  fullFrac = 1.0 - mean([binEdgesMS(b + 1), minMS])  * (1.0 - preStimMaxFrac) / (maxMS - minMS); % frac full trials in changing part of bin
  adjust = adjust + fullFrac * (binEdgesMS(b + 1) - minMS) / diff(binEdgesMS(b:b + 1));
  y(b) = y(b) / adjust;
  b = b + 1;
  
  % The remaining bins have a linearly changing number of short trials across their full width.
  for b = b:length(binEdgesMS) - 1
    fullFrac = 1.0 - (mean(binEdgesMS(b:b + 1)) - minMS) * (1.0 - preStimMaxFrac) / (maxMS - minMS); % average full trials across entire bin
    y(b) = y(b) / fullFrac;
  end
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

function doOnePlot(plotIndex, x, y, lineType, theTitle, theXLabel)

  subplot(4, 2, plotIndex);
  plot(x, y, lineType);
  title(theTitle);
  xlabel(theXLabel);
  xlim([0, inf]);
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
	headerText{length(headerText) + 1} = sprintf('RT window %d ms', file.rewardedLimitMS);
	headerText{length(headerText) + 1} = sprintf('Hit rate %.2f', pH);
  pFA = lambda * file.rewardedLimitMS;
	headerText{length(headerText) + 1} = sprintf('FA rate %.2f', pFA);
  [dP, c] = dPrime(pH, pFA);
	headerText{length(headerText) + 1} = sprintf('d-prime = %.2f; criterion %.2f', dP, c);

  text(0.00, 1.00, headerText, 'VerticalAlignment', 'top');
end