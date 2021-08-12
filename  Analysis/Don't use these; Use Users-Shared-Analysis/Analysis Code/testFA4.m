function testFA4
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
  rateFAs = [0.010, 0.05, 0.1, 0.5, 1.0]; % rate of 1/s for a time period of 1 s
  pHits = 0.50;
  file.preStimMaxMS = 3000;                % maximum prestim time
  file.rewardedLimitMS = 700;             % reaction time window
  repsPerLambda = 25;                     % number of repetitions
  numTrials = 500;
  figure(1);
  set(gcf, 'units', 'inches', 'position', [27, 10.0, 7.5, 10]);
  clf;  
  set(gca,'XScale','log','YScale','log');
  hold on;
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
  hold on;
  plot([rateFAs(1), rateFAs(end)], [rateFAs(1), rateFAs(end)], ':k');
  text(0.05, 0.95, {sprintf('maxStimMS %d', file.preStimMaxMS),...
    sprintf('%d Trials', numTrials), sprintf('%d fits per Lambda', repsPerLambda)}, 'Units', 'Normalized', ...
    'verticalAlignment', 'top', 'fontSize', 14);
  xlabel('Actual Rates');
  ylabel('Fitted Rates');
end

%%
function [rateS] = earlyRate(file, trials, indices)
  
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
  pastPreStim = FTrialTimesMS >= 0 & FTrialTimesMS < 3000;
  fTrialsTimesMS = FTrialTimesMS(pastPreStim);
  
  sumTimeMS = sum(fTrialsTimesMS) + sum([stimOnHMS, stimOnMMS]);
  L = (sumTimeMS / sum(pastPreStim) / 1000.0);
%   jS = file.preStimMaxMS / 1000.0;
%  
%   lambda = L / (1.0 - exp(-jS / L));
  fprintf('L rate %.4f ', 1.0 / L);
  rateS = 1 / L;
end

%%
function trials = makeTrials(rateS, pHit, numTrials, file)
%
% Make fake trials with given FA and H rates, with stim On uniformly between preStimMinMS and preStimMaxMS
%
  pF = 1.0 - exp(-rateS / 1000);          % Poisson probability for a false alarm (per millisecond)
  fCount = 0;
  for n = numTrials:-1:1                  % backward, to preallocate structure array
    trials(n).trialEnd = [];
    trial.preStimMS = (n - 1) * file.preStimMaxMS / (numTrials - 1);
    %{
    Step through time assigning FAs at random.  Note that the start time is important for matching event rates.  
    The events should start precisely at the time when counting is going to start, not from before or after.  
    This is not a factor when processing real data. 
    %}
    for t = 0:trial.preStimMS  
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
        trials(n).reactTimeMS = trial.preStimMS + 100;
      else
        trials(n).trialEnd = 2;             % miss trial
        trials(n).reactTimeMS = trial.preStimMS + file.rewardedLimitMS + 100;
      end
    end
    trials(n).trial = trial;                     % allocate in reverse to preallocate
  end
  fprintf('%d FAs\n', fCount);
end