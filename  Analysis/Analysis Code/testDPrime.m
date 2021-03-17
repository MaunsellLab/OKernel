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
%   OKMatlab([], file, trials);                             % display the results
  
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
  MeanPFAPerInc = mean(n ./ FANumTrials);
  pFAPerMS = 1.0 - exp(log((1.0 - MeanPFAPerInc)) / tInc);
  pFA = 1.0 - (1.0 - pFAPerMS)^(file.rewardedLimitMS - tooFastMS);
  
  % compile the probability of a false alarm for stimOns at different trial times
  fracFA = zeros(1, file.preStimMaxMS);
  trialTimePFA = zeros(1, file.preStimMaxMS);
  for t = 1:file.preStimMaxMS
    fracFA(t) = sum((FATrialTimeMS + tooFastMS) > t & FATrialTimeMS <= t + file.rewardedLimitMS) / numFAs;
    trialTimePFA(t) = sum((FATrialTimeMS + tooFastMS) > t & FATrialTimeMS <= t + file.rewardedLimitMS) / ...
      sum((FATrialTimeMS + tooFastMS) > t);
  end
  
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
  
  figure(1);
  clf;
    
  doOnePlot(1, edges, pHit, 'o', 'Hits/(Hits+Misses)', 'Time After StimOn');
  doOnePlot(4, 1:length(trialTimePFA), trialTimePFA, '-', 'FA Conditional Prob. RT Window Hit', 'Time After Trial Start');
  
  doOneHist(1, stimOnHMS, sprintf('Hits (n=%d)', length(indices.correct)), 'StimOn Time');
  doOneHist(2, stimOnMMS, sprintf('Misses (n=%d)', length(indices.miss)), 'StimOn Time');
  doOneHist(3, stimOnFAMS, sprintf('FAs (n=%d)', numFAs), '(Unreached) StimOn for FA');
  doOneHist(4, FATrialTimeMS, sprintf('FAs (n=%d)', numFAs), 'Lever Release for FA');
end

function doOneHist(row, values, theTitle, theXLabel)

  subplot(4, 2, (row - 1) * 2 + 1);
  histogram(values, 'binwidth', 250);
  title(theTitle);
  xlabel(theXLabel);
  xlim([0, inf]);
end

function doOnePlot(row, x, y, lineType, theTitle, theXLabel)

  subplot(4, 2, row * 2);
  plot(x, y, lineType);
  title(theTitle);
  xlabel(theXLabel);
  xlim([0, inf]);
  ylim([0, 1.0]);
end