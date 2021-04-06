%%
function earlyRate = earlyRate(file, trials, correctIndices, failIndices, earlyIndices)
%
% Calculate the average rate at which false alarms occurred within a data file.  The data must be from a task
% in which stimulus on times were uniformly distributed within trials.  It is assumed that false alarms occur
% at a constant rate, following Poisson statistics.
%
  % If no arguments are supplied, run a test
  if nargin == 0
    testEarlyRate();
    return;
  end
  
	% Extract trial outcomes.
  trialStructs = [trials(:).trial];                       % trial structs extracted from trials array
  preStimMS = [trialStructs(:).preStimMS];                % planned stimulus on time
  stimOnCMS = preStimMS(correctIndices);                  % stimOn for hits
  stimOnFMS = preStimMS(failIndices);
  stimOnEMS = preStimMS(earlyIndices);
  
  % Get the early times in trial using the RT and (tentative) preStimMS.
  earlyRTs = oneRTPerTrial(trials, earlyIndices);
  earlyTrialMS = stimOnEMS + earlyRTs;                    % time in trial when an early occurred
  leadS = file.preStimMinMS + file.tooFastMS;             % no responses count before the earliest stimOn + tooFast
  earlyPastMinTime = earlyTrialMS >= leadS;
  earlyTrialMS = earlyTrialMS(earlyPastMinTime) - leadS;
  stimOnCMS = stimOnCMS(stimOnCMS >= leadS) - file.preStimMinMS;
  stimOnFMS = stimOnFMS(stimOnFMS >= leadS) - file.preStimMinMS;
  sumTimeMS = sum([earlyTrialMS, stimOnCMS, stimOnFMS]);
  earlyRate = 1000.0 / (sumTimeMS / sum(earlyPastMinTime));
end

%%
function testEarlyRate
%{
Simulate Poisson false alarms and reconstruct the rate parameter.  Poisson events have delays that are distributed
exponentially, so fitting an exponential to the distribution yields the rate parameter.  Many of our tasks involve
stimuli appearing after uniformly-distributed random delays.  False alarms can only occur before those stimuli, so
stimulus occurrences effectively truncate the detection of later false alarms.  Because the random stimulus timing
is is uniformly distributed, it is straightforward to compensate for these truncated trials.
%}
  earlyRates = [0.003, 0.010, 0.03, 0.1, 0.3, 1.0, 3, 10]; % rate of 1/s for a time period of 1 s
  pHits = 0.5;
  file.preStimMinMS = 500;                % minimum prestim time
  file.preStimMaxMS = 3000;               % maximum prestim time
  file.tooFastMS = 100;                   % too fast time
  file.rewardedLimitMS = 700;             % reaction time window
  repsPerLambda = 25;                     % number of repetitions
  numTrials = 500;

  figure(10);
  set(gcf, 'units', 'inches', 'position', [27, 10.0, 7.5, 10]);
  clf;
  set(gca,'XScale','log','YScale','log');
  hold on;
  
  for p = 1:length(pHits)
    percentiles = zeros(3, length(earlyRates));
    for r = 1:length(earlyRates)
      estimates = zeros(1, repsPerLambda);
      for i = 1:repsPerLambda
        trials = makeTrials(earlyRates(r), pHits(p), numTrials, file);
        coreectIndices = find([trials(:).trialEnd] == 0);           % correct trials
        earlyIndices = find([trials(:).trialEnd] == 1);            % false alarm trials
        failIndices = find([trials(:).trialEnd] == 2);          % miss trials
        estimates(i) = earlyRate(file, trials, coreectIndices, failIndices, earlyIndices);
      end
      percentiles(:, r) = prctile(estimates, [25, 50, 75]);
    end
    errorbar(earlyRates, percentiles(2, :), percentiles(2, :) - percentiles(1, :), percentiles(2, :) ...
        - percentiles(3, :), 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k');
    hold on;
  end
  plot([earlyRates(1), earlyRates(end)], [earlyRates(1), earlyRates(end)], ':k');
  text(0.05, 0.95, {sprintf('minStimMS %d', file.preStimMinMS), sprintf('maxStimMS %d', file.preStimMaxMS),...
    sprintf('tooFastMS %d', file.tooFastMS) ...
    sprintf('%d Trials', numTrials), sprintf('%d fits per rate', repsPerLambda)}, 'Units', 'Normalized', ...
    'verticalAlignment', 'top', 'fontSize', 14);
  xlabel('Actual Early Rates');
  ylabel('Fitted Early Rates');
end

%%
function trials = makeTrials(rateS, pHit, numTrials, file)
%
% Make fake trials with given FA and H rates, with stimOn uniformly distributed between preStimMinMS and preStimMaxMS
%
  pEarly = 1.0 - exp(-rateS / 1000);          % Poisson probability for a false alarm (per millisecond)
  for n = numTrials:-1:1                  % backward, to preallocate structure array
    trials(n).trialEnd = [];
    trial.preStimMS = file.preStimMinMS + (n - 1) * (file.preStimMaxMS - file.preStimMinMS) / (numTrials - 1);
    %{
    Step through time assigning earlies at random.  Note that the start time is important for matching event rates.  
    The events should start precisely at the time when counting is going to start, not from before or after.  
    This is not a factor when processing real data. 
    %}
    for t = file.preStimMinMS + file.tooFastMS:trial.preStimMS + file.tooFastMS  
      if rand() < pEarly
        trials(n).trialEnd = 1;             % fa trial
        trials(n).reactTimeMS = t - trial.preStimMS;
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
    trials(n).trial = trial;               	% allocate in reverse to preallocate
  end
end