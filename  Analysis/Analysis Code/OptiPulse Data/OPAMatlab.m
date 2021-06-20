function OPAMatlab(file, trials, plotIndex)

% OPAMatlab is invoked at the end of every trial.  As a free-standing function that is instantiated on each call,
% it has no way to store static across calls.  Instead, such values are stored in a struct, dParams.  dParams
% arrives as an empty matrix on the first call, so the first call can be identified in this way.  By returning
% dParam with essential values, they can be recovered in each call.

  if nargin == 0
    OPExamples();
    return;
  end
  if ~isscalar(file.subjectNumber)
    file.subjectNumber = file.subjectNumber(end);
  end
  file.trials = size(trials, 2);              % update the number of trials
  if length([trials(:).trialEnd]) ~= length(trials)
    fprintf('Adding trialEnd to trial(s) without');
    for i = 1:length(trials)
      if isempty(trials(i).trialEnd)
          trials(i).trialEnd = -1;
      end
    end
  end

  trialStructs = [trials(:).trial];                  	% trial structs extracted from trials array
  psychometric(trials, trialStructs, plotIndex);
  row = struct;
  [row, indices, ~, stimIndices, noStimIndices] = getDPrimes(file, trials, row);
  drawText(file, trials, indices, stimIndices, noStimIndices, row, plotIndex);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performance Text

function drawText(file, trials, indices, stimIndices, noStimIndices, row, plotIndex)

  axisHandle = subplot(4, 3, plotIndex);						% default axes are 0 to 1
  hold off;
  set(axisHandle, 'visible', 'off');
  if plotIndex == 1
    text(0.00, 1.35, 'OptiPulse2', 'FontWeight', 'bold', 'FontSize', 16);
    text(0.00, 1.21, sprintf('Subject: %d; %s', file.subjectNumber, datestr(file.startTimeVec,'mmmm dd, yyyy HH:MM')), ...
      'FontSize', 12);
  end
  text(0.00, 1.10, datestr(file.startTimeVec,'mmmm dd, yyyy HH:MM'));
  text(0.00, 1.00, {' ', 'Trials:', 'Correct:', 'Failed:', 'Early:', 'd'''}, ...
    'verticalAlignment', 'top', 'horizontalAlignment', 'left');
  oneColumn(indices, noStimIndices | stimIndices, 'All', row.dPrime, 1);
  oneColumn(indices, noStimIndices, 'NoStim', row.noStimDPrime, 2);
  oneColumn(indices, stimIndices, 'Stim', row.stimDPrime, 3);
  set(gcf, 'visible', 'on');
  i = find([trials(:).powerMW] > 0, 1);
  text(-0.02, 0.34, {sprintf('delta d-Prime: %4.2f', row.stimDPrime - row.noStimDPrime) ...
    sprintf('power (mW): %.2f', trials(i).powerMW), ...
    sprintf('pulse dur (ms): %.0f', trials(i).trial.pulseDurMS), ...
    sprintf('pulse delay (ms): %.0f', trials(i).trial.delayMS)}, ...
    'verticalAlignment', 'top', 'horizontalAlignment', 'left');
end

function textHandle = oneColumn(indices, select, headerText, dPrime, offset)
	numCorrect = sum(indices.correct(select));
  numFailed = sum(indices.fail(select));
  numEarly = sum(indices.early(select));
  totalTrials = numCorrect + numFailed + numEarly;
  textHandle = text(0.30 * offset, 1.00, {headerText, sprintf('%.0f', totalTrials'), ...
        sprintf('%.0f%%', numCorrect / totalTrials * 100.0), ...
        sprintf('%.0f%%', numFailed / totalTrials * 100.0), ...
        sprintf('%.0f%%', numEarly / totalTrials * 100.0), ...
        sprintf('%.2f', dPrime) ...
        }); 
  set(textHandle, 'verticalAlignment', 'top', 'horizontalAlignment', 'left');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Psychometric Function %%%

function psychometric(trials, trialStructs, plotIndex)

  hit = [trials(:).trialEnd] == 0;
  miss = [trials(:).trialEnd] == 2;
  stimValueSet = unique([trials(hit | miss).visualStimValue]);
  numStim = length(stimValueSet);
  delay0 = [trialStructs(:).delayIndex] == 1;
  delay0Hits = zeros(1, numStim);
  delay0N = zeros(1, numStim);
  delay1 = [trialStructs(:).delayIndex] == 2;
  delay1Hits = zeros(1, numStim);
  delay1N = zeros(1, numStim);
  noStimHits = zeros(1, numStim);
  noStimN = zeros(1, numStim);
  for s = 1:numStim                                          % for each stim value
      stimTrial = [trials(:).visualStimValue] == stimValueSet(s);
      delay0Hits(s) = sum(stimTrial & hit & delay0);
      delay0N(s) = delay0Hits(s) + sum(stimTrial & miss & delay0);
      delay1Hits(s) = sum(stimTrial & hit & delay1);
      delay1N(s) = delay1Hits(s) + sum(stimTrial & miss & delay1);
      noStimHits(s) = sum(stimTrial & hit & ~delay0 & ~delay1);
      noStimN(s) = noStimHits(s) + sum(stimTrial & miss & ~delay0 & ~delay1);
  end
  [delay0HitRate, delay0pci] = binofit(delay0Hits, delay0N);
  [delay1HitRate, delay1pci] = binofit(delay1Hits, delay1N);
  [noStimHitRate, noStimPci] = binofit(noStimHits, noStimN);
  delay0YNeg = delay0HitRate - delay0pci(:, 1)';
  delay1YNeg = delay1HitRate - delay1pci(:, 1)';
  noStimYNeg = noStimHitRate - noStimPci(:, 1)';
  delay0YPos = delay0pci(:, 2)' - delay0HitRate;
  delay1YPos = delay1pci(:, 2)' - delay1HitRate;
  noStimYPos = noStimPci(:, 2)' - noStimHitRate;

  subplot(4, 3, plotIndex + 3);
  cla;
  hold off;
  errorbar(stimValueSet, noStimHitRate, noStimYNeg, noStimYPos, '-s', 'markersize', 6, ...
            'color', 'k', 'markerfacecolor', [0.75, 0.75, 0.75]);
  hold on;
  errorbar(stimValueSet, delay0HitRate, delay0YNeg, delay0YPos, '-s', 'markersize', 6, ...
           'color', [0.250, 0.250, 1.0], 'markerfacecolor', [0.50, 0.50, 1.0]);
  % hold on;
  errorbar(stimValueSet, delay1HitRate, delay1YNeg, delay1YPos, '-s', 'markersize', 6, ...
           'color', [0.250, 0.250, 1.0], 'markerfacecolor', [0.50, 0.50, 1.0]);
  set(gca,'xscale','log', 'xgrid', 'on');
  hold on;
  set(gca, 'ylim', [0 1]);
  xLimits = get(gca, 'XLim');
  set(gca, 'XLim', [5, min(xLimits(2), 100)]);
  set(gca, 'XTickLabel',num2str(transpose(get(gca, 'XTick'))))          % get rid of scientific notation
  ylabel('Percent Correct');
  xlabel('Gabor Contrast (%)');
end
 
 %%
function [row, theIndices, trials, stimIndices, noStimIndices] = getDPrimes(file, trials, row)
%
% Compute d-Primes for one session
%
% To get the d-Primes, we need to find the stimulus contrast that was stimulated, and then we need to compare it with
% only unstimulated trials that used the same visual contrast.

  powerMW = [trials(:).powerMW];                        % get power applied for each trial                        
  stimIndices = powerMW > 0;
  if sum(stimIndices) == 0                                    % no trials with full opto contrast
    error('No stimulated files found');
  end
  visualStimValues = [trials(:).visualStimValue];
  contrastPC = unique(visualStimValues(stimIndices));
  if length(contrastPC) > 1
    error('Multiple contrasts stimulated');
  end
  noStimIndices = ~stimIndices & visualStimValues == contrastPC;
  
  % only consider the range between the first and last stimulated trials -- trim the vectors
  firstStimIndex = find(stimIndices > 0, 1);                 	% first stimulated trial
	lastStimIndex = find(stimIndices > 0, 1, 'last');          	% last stimulated trial
  stimIndices = stimIndices(firstStimIndex:lastStimIndex);
  noStimIndices = noStimIndices(firstStimIndex:lastStimIndex);
  trials = trials(firstStimIndex:lastStimIndex);

  % consider only the trials that used the stimulated visual stimulus contrast
  usableIndices = stimIndices | noStimIndices;
  trials = trials(usableIndices);
  stimIndices = stimIndices(usableIndices);
  noStimIndices = noStimIndices(usableIndices);
  eotCodes = zeros(1, length(trials));                        % one and only one RT and endTrial per trial
  for t = 1:length(trials)                               
    if length(trials(t).reactTimeMS) > 1
        trials(t).reactTimeMS = trials(t).reactTimeMS(1);   
    end
    if ~isfield(trials(t), 'trialEnd') || isempty(trials(t).trialEnd)
        trials(t).trialEnd = -1;
    end
    if length(trials(t).trialEnd) > 1
        trials(t).trialEnd = trials(t).trialEnd(1);   
    end
    eotCodes(t) = trials(t).trialEnd;
  end

  % calculate the overall d'
  theIndices = allIndices(trials, eotCodes);
  if sum(theIndices.correct | theIndices.fail | theIndices.early) < 10
    row = [];
    return;
  end
  % find the response interval and get a modified set of indices that limits hits to only that interval
	[respLimitsMS, theIndices, ~, ~] = getResponseLimits(file, trials, theIndices);
  row.startRT = respLimitsMS(1);
  row.endRT = respLimitsMS(2);
  row.RTWindowMS = diff(respLimitsMS);
  row.noStimCorrects = sum(theIndices.correct & ~stimIndices);
  row.noStimFails = sum(theIndices.fail & ~stimIndices);
  row.noStimEarlies = sum(theIndices.early & ~stimIndices); 
  row.numNoStim = row.noStimCorrects + row.noStimFails + row.noStimEarlies;  % don't count kEOTIgnored
  row.stimCorrects = sum(theIndices.correct & stimIndices);
  row.stimFails = sum(theIndices.fail & stimIndices);
  row.stimEarlies = sum(theIndices.early & stimIndices);
  row.numStim = row.stimCorrects + row.stimFails + row.stimEarlies;                  % don't count kEOTIgnored

  % find the performance across stim and nostim trials combined
  hitRate = sum(theIndices.correct) / (sum(theIndices.correct) + sum(theIndices.fail));
  rateEarly = earlyRate(file, trials, theIndices.correct, theIndices.fail, theIndices.early);
  row.pFA = 1.0 - exp(-rateEarly * row.RTWindowMS / 1000.0);
  row.pHit = (hitRate - row.pFA) / (1.0 - row.pFA);
  [row.dPrime, row.c] = dprime(row.pHit, row.pFA, true);
  if abs(row.dPrime) == Inf
    row.dPrime = NaN; 
  end

  % calculate the nostim trial d' using the all-trial pFA and the all-trial response window
  indices.correct = theIndices.correct & ~stimIndices;
  indices.fail = theIndices.fail & ~stimIndices;
  indices.early = theIndices.early & ~stimIndices;
  if sum(indices.correct | indices.fail | indices.early) > 0
    hitRate = sum(indices.correct) / (sum(indices.correct) + sum(indices.fail));
    rateEarly = earlyRate(file, trials, indices.correct, indices.fail, indices.early);
    row.noStimPFA = 1.0 - exp(-rateEarly * row.RTWindowMS / 1000.0);
    row.noStimPHit = (hitRate - row.pFA) / (1.0 - row.pFA);         % using overall pFA
    [row.noStimDPrime, row.noStimC] = dprime(row.noStimPHit, row.pFA, true);
    if abs(row.noStimDPrime) == Inf
      row.noStimDPrime = NaN; 
    end
  end
    
  % calculate the stim trial d' using the all-trial pFA and the all-trial response window
  indices.correct = theIndices.correct & stimIndices;
  indices.fail = theIndices.fail & stimIndices;
  indices.early = theIndices.early & stimIndices;
  if sum(indices.correct | indices.fail | indices.early) > 0
    hitRate = sum(indices.correct) / (sum(indices.correct) + sum(indices.fail));
    rateEarly = earlyRate(file, trials, indices.correct, indices.fail, indices.early);
    row.stimPFA = 1.0 - exp(-rateEarly * row.RTWindowMS / 1000.0);
    row.stimPHit = (hitRate - row.pFA) / (1.0 - row.pFA);   % using overall pFA
    [row.stimDPrime, row.stimC] = dprime(row.stimPHit, row.pFA, true);
    if abs(row.stimDPrime) == Inf
      row.stimDPrime = NaN; 
    end
  end
end