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

  load('/Users/Shared/Data/OKernel/902/MatFiles/2019-09-16.mat', 'file', 'trials');
  OKMatlab([], file, trials);                             % display the results
  
  indices.correct = find([trials(:).trialEnd] == 0);      % correct trials
  indices.fa = find([trials(:).trialEnd] == 1);           % false alarm trials
  indices.miss = find([trials(:).trialEnd] == 2);         % miss trials
  numFAs = length(indices.fa);
  trialStructs = [trials(:).trial];                       % trial structs extracted from trials array
  preStimMS = [trialStructs(:).preStimMS];                % planned stimulus on time
  tooFastMS = file.tooFastMS;
  
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
  
% compile the probability of a false alarm for stimOns at different trial times
%   MeanPFAPerInc = mean(n ./ FANumTrials);
%   pFAPerMS = 1.0 - exp(log((1.0 - MeanPFAPerInc)) / tInc);
%   pFA = 1.0 - (1.0 - pFAPerMS)^(file.rewardedLimitMS - tooFastMS);
%   fracFA = zeros(1, file.preStimMaxMS);
%   trialTimePFA = zeros(1, file.preStimMaxMS);
%   for t = 1:file.preStimMaxMS
%     fracFA(t) = sum((FATrialTimeMS + tooFastMS) > t & FATrialTimeMS <= t + file.rewardedLimitMS) / numFAs;
%     trialTimePFA(t) = sum((FATrialTimeMS + tooFastMS) > t & FATrialTimeMS <= t + file.rewardedLimitMS) / ...
%       sum((FATrialTimeMS + tooFastMS) > t);
%   end
  
  % compile the hit probability for stimOns in different bins during the stimOn interval
  hits = zeros(1, length(edges));
  misses = zeros(1, length(edges));
  FAs = zeros(1, length(edges));
  CRs = zeros(1, length(edges));
%   for h = 1:length(indices.correct)
%     stimOnTimeMS = preStimMS(indices.correct(h));
%     edgeBin = floor(stimOnTimeMS / tInc);
%     hits(edgeBin) = hits(edgeBin) + 1;              % add one hit trial to this bin
%     numFACR = sum(FATrialTimeMS > stimOnTimeMS + tooFastMS); % number of non H/M trials in play for stimOn time
%     fFA = fracFA(stimOnTimeMS);                     % fraction of FAs that would have hit this RT window
%     nFA = fFA * numFAs;                             % number of FA trials that would have hit this RT window
%     nCR = numFAs - nFA;                             % number of CR trials 
%   end
  
  % compile a histogram of the hit rate as a function of stimOn time.
  pHit = zeros(1, length(edges));
  for t = 1:length(edges)
    h = sum(stimOnHMS >= edges(t) & stimOnHMS < edges(t) + tInc);
    m = sum(stimOnMMS >= edges(t) & stimOnMMS < edges(t) + tInc);
    pHit(t) = h / (h + m);
  end
  
%   dprime(mean(pHit), pFA)
  
  f2 = figure(2);
  clf;
  set(gcf, 'units', 'inches', 'position', [27, 10.0, 7.5, 10]);    
  doOnePlot(2, edges, pHit, 'o', 'Hits/(Hits+Misses)', 'Time After StimOn');
%   doOnePlot(4, 1:length(trialTimePFA), trialTimePFA, '-', 'FA Conditional Prob. RT Window Hit', 'Time After Trial Start');
  
  doOneHist(1, stimOnHMS, sprintf('Hits (n=%d)', length(indices.correct)), 'StimOn Time');
  doOneHist(3, stimOnMMS, sprintf('Misses (n=%d)', length(indices.miss)), 'StimOn Time');
  doOneHist(5, stimOnFAMS, sprintf('FAs (n=%d)', numFAs), '(Unreached) StimOn for FA');
%   doOneHist(7, FATrialTimeMS, sprintf('FAs (n=%d)', numFAs), 'Lever Release for FA');
% 	doOneHist(6, [stimOnHMS, stimOnMMS, stimOnFAMS], sprintf('H+M+FA (n=%d)', ...
%     length(indices.correct) + length(indices.miss) + numFAs), 'StimOn Time (reached or unreached)');
	doOneHist(4, [stimOnHMS, stimOnMMS], sprintf('H+M (n=%d)', ...
    length(indices.correct) + length(indices.miss)), 'StimOn Time');

  
	numF = length(FATrialTimeMS);
  dT = 500;
 
  % Produce a linear fit to trials unavailable for false alarms because they were consumed by hits or misses. 
  % The normal assumption is that available trials go to zero at preStimMaxMS, but that is not always the case.  
  binEdgesMS = file.preStimMinMS:dT:file.preStimMaxMS;
	x = binEdgesMS(1:end - 1) + dT / 2.0;
  y = histcounts([stimOnHMS, stimOnMMS], binEdgesMS);   % combined hit and miss stimOn distributions
  coeff = polyfit(x, y, 1);
  preStimMaxFrac = (file.preStimMaxMS * coeff(1) + coeff(2)) / coeff(2);
  
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
  subplot(4, 2, 7);
  hold on;
  plot(fitted, x, y);
  hl = findobj(gcf,'type','legend');
  delete(hl);
  xlabel('Time (ms)');
  ylabel('');
  axis([0, file.preStimMaxMS, 0, inf]);
  text(0.95, 0.95, sprintf('L = %.5f', -fitted.b), 'units', 'normalized', 'horizontalAlignment', 'right',...
    'verticalAlignment', 'top', 'fontSize', 12);
  doOneHist(7, FATrialTimeMS, sprintf('FAs (n=%d)', numFAs), 'FA Time from Trial Start (ms)');

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
function doOneHist(plotIndex, values, theTitle, theXLabel)

  subplot(4, 2, plotIndex);
  histogram(values, 'binwidth', 250);
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