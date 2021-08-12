function testFA2
%{
Simulate Poisson false alarms and reconstruct the rate parameter.  Poisson events have delays that are distributed
exponentially, so fitting an exponential to the distribution yields the rate parameter.  Many of our tasks involve
stimuli appearing after uniformly-distributed random delays.  False alarms can only occur before those stimuli, so
stimulus occurrences effectively truncate the detection of later false alarms.  Because the random stimulus timing
is is uniformly distributed, it is straightforward to compensate for these truncated trials.

This test displays three runs of the simulation.  The first has full trials with no truncations.  The performance
is good.  The second has truncated, uncorrected trials, for which the rate is systematically overestimated, particularly
for low rates. The third has truncated trials with false alarm counts that have been compensated for the truncation.
%}
  rateFAs = [0.10, 0.25, 0.5, 0.75, 1.0];          % rate of 1/s for a time period of 1 s
  pHits = [0.0, 0.50, 1.0];
  file.preStimMinMS = 500;          % minimum prestim time
  file.preStimMaxMS = 3000;         % maximum prestim time
  file.tooFastMS = 100;             % too fast time
  file.rewardedLimitMS = 700;       % reaction time window
  repsPerLambda = 25;              % number of repetitions
  numTrials = 500;
  f1 = figure(1);
  set(gcf, 'units', 'inches', 'position', [27, 10.0, 7.5, 10]);
  clf;
  
  for p = 1:length(pHits)
    percentiles = zeros(3, length(rateFAs));
    for r = 1:length(rateFAs)
      estimates = zeros(1, repsPerLambda);
      for i = 1:repsPerLambda
        trials = makeTrials(rateFAs(r), pHits(p), numTrials, file);
        indices.correct = find([trials(:).trialEnd] == 0);      % correct trials
        indices.fa = find([trials(:).trialEnd] == 1);           % false alarm trials
        indices.miss = find([trials(:).trialEnd] == 2);         % miss trials
        estimates(i) = earlyRate(file, trials, indices);
      end
      percentiles(:, r) = prctile(estimates, [25, 50, 75]);
    end
    errorbar(rateFAs, percentiles(2, :), percentiles(2, :) - percentiles(1, :), percentiles(2, :) - percentiles(3, :), ...
      'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k');
    hold on;
  end
  a = axis;
  aMax = max(a(3:4));
  axis([0, aMax, 0, aMax]);
  hold on;
  plot([0, aMax], [0, aMax], ':k');
  text(0.05, 0.95, {sprintf('minStimMS %d', file.preStimMinMS), sprintf('maxStimMS %d', file.preStimMaxMS),...
    sprintf('%d Trials', numTrials), sprintf('%d fits per Lambda', repsPerLambda)}, 'Units', 'Normalized', ...
    'verticalAlignment', 'top', 'fontSize', 14);
  xlabel('Actual Rates');
  ylabel('Fitted Rates');
end

%%
function rateS = earlyRate(file, trials, indices)
  
  trialStructs = [trials(:).trial];                       % trial structs extracted from trials array
  preStimMS = [trialStructs(:).preStimMS];                % planned stimulus on time
  stimOnHMS = preStimMS(indices.correct);                 % stimOn for hits
  stimOnMMS = preStimMS(indices.miss);
  stimOnFMS = preStimMS(indices.fa);
  
  % Get the FA times in trial using the RT and (tentative) preStimMS, making sure that there is only one RT per trial.
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
  
% The times after hit/miss truncations need to be incremented in proportion to the fraction of trials
% that have been truncated by hit/miss trials at that point. The fractional increment is based on a count
% that excludes all trial for which a FA has already occurred for abbreviated trials
%
  leadS = file.preStimMinMS + file.tooFastMS;
  pastPreStim = FTrialTimesMS >= leadS;
  sumTimeMS = sum(FTrialTimesMS(pastPreStim) - leadS) + sum([stimOnHMS, stimOnMMS] - leadS);
  lambdaS = sumTimeMS / sum(pastPreStim) / 1000.0;
  
%{
  We need to correct for the fact that we didn't allow for intervals greater than preStimMaxMS + tooFastMS.
  If we flip the axes on the exponential that defines the distribution of time to FA (y = exp(-L * x), we get 
  y = -1/L * log(x) (where L is lambda). We can get the average value of y (the FA period) by evaluating the
  integral over the interval from 0 to 1 ( = -x/L * (log(x) - 1), which reduces to 1/L over the interval 0-1).
  
  This integral assumes that the distribution extends to y = inf, but in fact we truncated it at 
  y = preStimMaxMS + tooFastMS, so the estimate of L is too large. We can correct this by finding the truncated
  distribution that has our average value. The truncation occurs at y = file.preStimMaxMS + file.tooFastMS.  This
  corresponds to an x value of xT = exp(-L * (file.preStimMaxMS + file.tooFastMS) / 1000.0).  To get the area of the
  truncated distribution, we remove the area in the exponential distribution from 0 to xT, then add back a retangular
  area corresponding to the xT * (file.preStimMaxMS + file.tooFastMS).
  
  The areas of in the exponential distribution from 0 to xT is -xT/L * (log(xT) - 1).  So we have the full distribution
    1.0 / L, minus the truncated portion, -xT/L * (log(xT) - 1) plus the rectangular region. The first to reduce to
    1.0 / L * (1.0 - xT * (log(xT) - 1))
%}
  integral = 1.0 / lambdaS;
  truncX = exp(-lambdaS * (file.preStimMaxMS + file.tooFastMS) / 1000.0);
  truncIntegral = (1.0 / lambdaS) + truncX / lambdaS * (log(truncX) - 1) ...
        + truncX * (file.preStimMaxMS + file.tooFastMS) / 1000.0;
  while truncIntegral < integral
    lambdaS = lambdaS - 0.001;
    truncX = exp(-lambdaS * (file.preStimMaxMS + file.tooFastMS) / 1000.0);
    truncIntegral = 1.0 / lambdaS ...
        + truncX / lambdaS * (log(truncX) - truncX) ...
        + truncX * (file.preStimMaxMS + file.tooFastMS) / 1000.0;
  end
  rateS = 1.0 / lambdaS;
end

%%
function trials = makeTrials(rateS, pHit, numTrials, file)
%
% Make fake trials with given FA and H rates, with stim On uniformly between preStimMinMS and preStimMaxMS
%
  fCount = 0;
  pF = 1.0 - exp(-rateS / 1000);          % Poisson probability for a false alarm (per millisecond)
  for n = numTrials:-1:1                  % backward, to preallocate structure array
    trials(n).trialEnd = [];
    trial.preStimMS = file.preStimMinMS + (n - 1) * (file.preStimMaxMS - file.preStimMinMS) / (numTrials - 1);
    %{
    Step through time assigning FAs at random.  Note that the start time is important for matching event rates.  
    The events should start precisely at the time when counting is going to start, not from before or after.  
    This is not a factor when processing real data. 
    %}
    for t = file.preStimMinMS + file.tooFastMS:trial.preStimMS + file.tooFastMS  
      if rand() < pF
        trials(n).trialEnd = 1;             % fa trial
        trials(n).reactTimeMS = t - trial.preStimMS;
        fCount = fCount + 1;
        break;
      end
    end
    if isempty(trials(n).trialEnd) 
      if rand() < pHit
        trials(n).trialEnd = 0;             % hit trial
        trials(n).reactTimeMS = trial.preStimMS + file.tooFastMS + 100;
      else
        trials(n).trialEnd = 2;             % miss trial
        trials(n).reactTimeMS = trial.preStimMS + file.rewardedLimitMS + 100;
      end
    end
    trials(n).trial = trial;                     % allocate in reverse to preallocate
  end
end