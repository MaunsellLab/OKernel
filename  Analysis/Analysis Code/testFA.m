function testFA
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
  minMS = 500;          % minimum prestim time
  maxMS = 3000;         % maximum prestim time
  repsPerLambda = 100;  % number of repetitions
  numTrials = 5000;
  dT = 500;
  lambdas = [0.0005:0.00025:0.0020];      % rate of 1/s for a time period of 1 ms
  funcs = {@doFullFit, @doTruncFit, @doTruncCorrFit};
  titles = {'Full Length Trials', 'Trials Truncated by StimOn', ...
    'Trials Truncated by StimOn and Corrected for Truncation'};
  f1 = figure(1);
  set(gcf, 'units', 'inches', 'position', [27, 10.0, 7.5, 10]);
  clf;
  f2 = figure(2);
  set(gcf, 'Units', 'inches', 'Position', [27, 0.5, 7.5, 10]);
  clf;
  
  numFunctions = length(funcs);
  meanFits = zeros(numFunctions, length(lambdas));
  medianFits = zeros(numFunctions, length(lambdas));
  quartiles = zeros(numFunctions, length(lambdas), 2);
  for f = 1:numFunctions
    for l = 1:length(lambdas)
      fitL = zeros(1, repsPerLambda);
      numF = zeros(1, repsPerLambda);
      for r = 1:repsPerLambda
        while numF(r) == 0                  % zero returned when there are too few points to fit
          [fitted, numF(r), x, y] = funcs{f}(lambdas(l), numTrials, minMS, maxMS, dT);
        end
        fitL(r) = -fitted.b;
        if f == 1 && r == 1
          figure(f2);
          subplot(4, 3, l);
          plot(fitted, x, y);
          hl = findobj(gcf,'type','legend');
          delete(hl);
          xlabel('Time (ms)');
          ylabel('');
          xlim([0, maxMS]);
          text(0.95, 0.95, sprintf('L = %.5f', lambdas(l)), 'Units', 'Normalized', 'horizontalAlignment', 'right',...
            'verticalAlignment', 'top', 'fontSize', 12);
        end
      end
      meanFits(f, l) = mean(fitL);
      medianFits(f, l) = median(fitL);
      quartiles(f, l, :) = prctile(fitL, [25, 75]);
      fprintf(' Actual %.5f, median (IQR) fit = %.5f (%.5f-%0.5f) (average Fs = %.1f)\n', ...
        lambdas(l), medianFits(f, l), quartiles(f, l,:), mean(numF));
    end

    figure(f1);
    subplot(numFunctions, 1, f);
    errorbar(lambdas, medianFits(f, :), medianFits(f, :) - quartiles(f, :, 1), quartiles(f, :, 2) - medianFits(f, :), ...
      'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k');
    a = axis;
    aMax = max(a(3:4));
    axis([0, aMax, 0, aMax]);
    hold on;
    h1 = plot([0, aMax], [0, aMax], ':k');
    text(0.05, 0.95, {sprintf('minStimMS %d', minMS), sprintf('maxStimMS %d', maxMS), sprintf('%d Trials', numTrials),...
      sprintf('%d fits per Lambda', repsPerLambda)}, 'Units', 'Normalized', 'verticalAlignment', 'top', ...
      'fontSize', 14);
    xlabel('Actual Lambda');
    ylabel('Fitted Lambda');
    title(titles{f});
  end
end

%%
function y = adjustForShortTrials(y, binEdgesMS, minMS, maxMS)
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
  adjust = (minMS - binEdgesMS(b)) / diff(binEdgesMS(b:b + 1));               % portion before minMS
  fullFrac = 1.0 - mean([binEdgesMS(b + 1), minMS]) / (maxMS - minMS);  % frac full trials in changing part of bin
  adjust = adjust + fullFrac * (binEdgesMS(b + 1) - minMS) / diff(binEdgesMS(b:b + 1));
  y(b) = y(b) / adjust;
  b = b + 1;
  
  % The remaining bins have a linearly changing number of short trials across their full width.
  for b = b:length(binEdgesMS) - 1
    fullFrac = 1.0 - (mean(binEdgesMS(b:b + 1)) - minMS) / (maxMS - minMS); % average full trials across entire bin
      y(b) = y(b) / fullFrac;
  end
end

%%
function [fitted, numF, x, y] = doFullFit(lambda, numTrials, ~, maxMS, dT)
%
% Do fits with trials that have no truncations
%
  numF = 0;
  pF = lambda * exp(-lambda); % Poisson probability for (1) F
  delaysMS = [];
  for n = 1:numTrials
    for t = 1:maxMS
      if rand() < pF
        delaysMS = [delaysMS, t]; %#ok<*AGROW>
        numF = numF + 1;
        break;
      end
    end
  end
  binEdgesMS = 0:dT:maxMS;
  y = histcounts(delaysMS, binEdgesMS); % counts adjusted for bin widths
  x = binEdgesMS(1:end - 1) + dT / 2.0;
  lastY = find(y > 0, 1, 'last' );
  y = y(1:lastY)';
  x = x(1:lastY)';
  fitted = fit(x, y, 'exp1', 'startPoint', [0, 0]);
end

%%
function [fitted, numF, x, y] = doTruncCorrFit(lambda, numTrials, minMS, maxMS, dT)
%
% Fits where the trials stop early, uniformly between minMS and maxMS
%
  pF = lambda * exp(-lambda); % Poisson probability for (1) F
  delaysMS = [];
  stimMSs = zeros(1, numTrials);
  for n = 1:numTrials
    stimMS = minMS + (n - 1) * (maxMS - minMS) / (numTrials - 1);
    stimMSs(n) = stimMS;                    % truncate trial at stimOn time
    for t = 1:stimMSs(n)
      if rand() < pF
        delaysMS = [delaysMS, t]; %#ok<*AGROW>
        break;
      end
    end
  end
  numF = length(delaysMS);
  binEdgesMS = 0:dT:maxMS;
  x = binEdgesMS(1:end - 1) + dT / 2.0;
  y = histcounts(delaysMS, binEdgesMS);
  
  % adjust for the abbreviated trials
  y = adjustForShortTrials(y, binEdgesMS, minMS, maxMS);
  
  % fit() doesn't work well if there are a lot of y == 0 entries at the end of the sequence.  These are undersampled,
  % and they can't be corrected for the number of short trials. 
  lastY = find(y > 0, 1, 'last' );
  y = y(1:lastY)';
  x = x(1:lastY)';
  fitted = fit(x, y, 'exp1', 'startPoint', [0, 0]);
end

%%
function [fitted, numF, x, y] = doTruncFit(lambda, numTrials, minMS, maxMS, dT)
%
% Fits where the trials stop early, uniformly between minMS and maxMS
%
  pF = lambda * exp(-lambda); % Poisson probability for (1) F
  delaysMS = [];
  stimMSs = zeros(1, numTrials);
  for n = 1:numTrials
    stimMS = minMS + (n - 1) * (maxMS - minMS) / (numTrials - 1);
    stimMSs(n) = stimMS;                    % truncate trial at stimOn time
    for t = 1:stimMSs(n)
      if rand() < pF
        delaysMS = [delaysMS, t]; %#ok<*AGROW>
        break;
      end
    end
  end
  numF = length(delaysMS);
  binEdgesMS = 0:dT:maxMS;
  x = binEdgesMS(1:end - 1) + dT / 2.0;
  y = histcounts(delaysMS, binEdgesMS);

  % fit() doesn't work well if there are a lot of y == 0 entries at the end of the sequence.  These are undersampled,
  % and they can't be corrected for the number of short trials. 
  lastY = find(y > 0, 1, 'last' );
  y = y(1:lastY)';
  x = x(1:lastY)';
  fitted = fit(x, y, 'exp1', 'startPoint', [0, 0]);
end