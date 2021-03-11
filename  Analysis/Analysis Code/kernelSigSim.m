function kernelSigSim()
% simulate the effects of selecting kernels with a prescribed level of
% significance

%{
This simulation tests whether our significance criterion creates kernels from noise.  The first figure simply plots
some representative kernels and histograms showing the distributions of the normalized power within them.  Gray bands
mark the 95% CIs on the kernels (binomial statistics), and dotted lines mark the same values on the histograms.  Red
lines on the histograms show the actual 5% and 95% quartiles for the distributions

The second plot shows the effects of selecting kernels.  The first distribution is the average of many session kernel
with no selection.  The second selects kernels that have a negative peak less than 1.25 * the CI in the 3rd quarter of
the kernel.  The third plot selects kernels that have a positive peak in the 3rd quarter.  The final plot selects 
kernels that have either a negative or positive peak in the 3rd quarter.  We select the same number of sessions for 
each test, so the CIs are all the same, and the expected variance is the same. Selecting positive or negative peaks
has the expected effect.  Selecting positive and negative peaks has no obvious effect.

The third plot checks whether the selection affects the variance of the computed kernels.  For each of the four
selection conditions, a variance is computed for the first and second half of the kernel.  This is based on the average
of many sessions, and the process is repeated many times.  For the unselected and pos/neg selection, there is no
obvious difference in the variance for the early (unselected) and late (selected portions).  For this analysis we
include the 3rd (selected) and 4th (unselected) quarters, to reveal the local varaince.  For this test, we do not use
simple variance because the offset means for the negative and positive selections will reduce the apparent variance. 
We instead keep the mean at zero and simply take the mean squared deviation from zero.
%}

  pHit = 0.50;      % hit rate
  nTrials = 200;    % number of trials per session
  nSessions = 100; % number of session
  nExamples = 8;    % number of examples to plot
  lenKernel = 800;  % length of kernel (nominal ms)
  titles = {'Unselected', '401-600 -1.5 CI', '401-600 +1.5 CI', '401-600 +/-1.5 CI'};
  conditions = {'normal', 'negative', 'positive', 'pos-neg'};

  % create a set of kernels
  [~, kernels] = makeKernels('normal', nSessions, nTrials, lenKernel, pHit);
  
  % plot some examples from the set of kernels
  % set up figure
  h = figure(1);
  set(h, 'Units', 'inches', 'Position', [25, 14.5, 10.24, 6.40]);
  clf;
  
  % plot the kernels and the CI
  CI = stimCI(nTrials);           % 95% CI for these kernels
  for s = 1:nExamples
    subplot(4, 4, s * 2 - 1);
    h = fill([0, size(kernels, 2), size(kernels, 2), 0], [CI, CI, -CI, -CI], [0.8, 0.8, 0.8]);
    set(h, 'linestyle', ':', 'facealpha', 0.25);
    hold on;
    plot(kernels(s, :), 'b');
    xlim([0, size(kernels, 2)]);
  end
  limits = sameYAxisScaling(4, 4, 1:2:nExamples * 2);
  
  % plot a histogram for each kernel, mark the CI
  limit = max(abs(limits));
  edges = linspace(-limit, limit, 12);
  q = zeros(nExamples, 2);
  for s = 1:nExamples
    subplot(4, 4, s * 2);
    h = histogram(kernels(s, :), edges);
    q(s, :) = quantile(h.Data, [0.05, 0.95]);
  end
  yLimits = sameYAxisScaling(4, 4, 2:2:nExamples * 2);
  for s = 1:nExamples
    subplot(4, 4, s * 2);
    hold on;
    plot([-CI, -CI], yLimits, 'k:');
    plot([CI, CI], yLimits, 'k:');
    plot([q(s, 1), q(s, 1)], yLimits, 'r:');
    plot([q(s, 2), q(s, 2)], yLimits, 'r:');
  end
  
  % set up figure for grand kernels
  h = figure(2);
  set(h, 'Units', 'inches', 'Position', [25, 7.5, 10.24, 6.40]);
  clf;
  
  % plot the grand kernel that includes all sessions
  CI = stimCI(nTrials * nSessions);           % 95% CI for grand kernel
  for c = 1:length(conditions)
    if c == 1
      titleStr = {sprintf('%d Sessions %d Trials', nSessions, nTrials), titles{1}};
    else
      titleStr = titles{c};
    end
    grandKernel = makeKernels(conditions{c}, nSessions, nTrials, lenKernel, pHit);
    plotGrandKernel(grandKernel, CI, c, titleStr, makeTextString(grandKernel));
  end
  sameYAxisScaling(2, 4, 1:4);
  
  % set up early/late variance scatterplots
  h = figure(3);
  set(h, 'Units', 'inches', 'Position', [25, 0.25, 10.24, 6.40]);
  clf;
  
  nReps = 250;
  for c = 1:length(conditions)
    if c == 2 || c == 3
      continue;
    end
    vars = zeros(nReps, 2);
    for r = 1:nReps
      fprintf('condition %d of %d; rep %d of %d\n', c, length(conditions), r, nReps);
      grandKernel = makeKernels(conditions{c}, nSessions, nTrials, lenKernel, pHit);
      vars(r, 1) = mean(sum(grandKernel(1:400).^2));
      vars(r, 2) = mean(sum(grandKernel(401:800).^2));
    end
    subplot(1, 4, c);
    plot(vars(:, 1), vars(:, 2), 'bo');
    hold on;
    axis square;
    yl = ylim();
    xl = xlim();
    maxVal = max(xl(2), yl(2));
    axis([0, maxVal, 0, maxVal]);
    plot([0, maxVal], [0, maxVal], 'k:');
    xlabel('Early');
    ylabel('Late');
    title(titles{c});
    text(0.05, 0.05, sprintf('E %.4f; L %.4f', mean(vars(:, 1)), mean(vars(:, 2))), 'units', 'normalized');
    drawnow;
  end
end

function [grandKernel, kernels] = makeKernels(condition, nSessions, nTrials, lenKernel, pHit)

  pulseWidth = 25;                          % width of pulse (nominal ms)
  testInterval = 401:600;                   % time interval to check for criterion
  kernels = zeros(nSessions, lenKernel);    % array of kernels to return
  critCI = 1.5 * stimCI(nTrials);           % 95% CI for these kernels
  
  for s = 1:nSessions
    while true
      bValues = randi([0, 1], nTrials, int32(lenKernel / pulseWidth));  % random binary values (one per pulse)
      optoStim = repelem(bValues, 1, pulseWidth);                      	% expand to fill pulsewidths

      % give each stimulus a random phase by circularly shifting by a random portion of pulse width
      for t = 1:nTrials
        optoStim(t, :) = circshift(optoStim(t, :), randi([0, pulseWidth - 1]), 2);
      end

      % create the kernel by splitting trials between hits and misses and subtracting the difference
      nHits = int32(nTrials * pHit);
      if pHit >= 1.0
        kernels(s, :) = sum(optostim) / nTrials;
      elseif pHit <= 0.0
        kernels(s, :) = -sum(optoStim) / nTrials;
      else
        kernels(s, :) = (sum(optoStim(1:nHits, :)) - sum(optoStim(nHits + 1:end, :))) / nTrials;
      end
      
      % check whether the current kernel satisfies the condition
      switch condition
        case 'normal'       % all kernels are accepted in the normal mode
          break;
        case 'negative'     % accepted only kernels that reach negative criterion during the test interval
          if min(kernels(s, testInterval), [], 2) < -critCI
            break;
          end
        case 'positive'     % accepted only kernels that reach negative criterion during the test interval
          if max(kernels(s, testInterval), [], 2) > critCI
            break;
          end
        case 'pos-neg'
          if min(kernels(s, testInterval), [], 2) < -critCI || max(kernels(s, testInterval), [], 2) > critCI
            break;
          end
        otherwise
          error('makeKernels -- unrecognized condition %s', condition);
      end
    end
  end
  grandKernel = mean(kernels);
end

function textStr = makeTextString(kernel)
  earlyVar = mean(sum(kernel(1:400).^2));
  lateVar = mean(sum(kernel(401:800).^2));
  textStr = {sprintf('early var: %.4f', earlyVar), sprintf(' late var: %.4f', lateVar)};
end

function plotGrandKernel(grandKernel, CI, plotIndex, titleStr, textStr)
  subplot(2, 4, plotIndex);
  h = fill([0, size(grandKernel, 2), size(grandKernel, 2), 0], [CI, CI, -CI, -CI], [0.8, 0.8, 0.8]);
  set(h, 'linestyle', ':', 'facealpha', 0.25);
  hold on;
  plot(grandKernel, 'b');
  xlim([0, size(grandKernel, 2)]);
  title(titleStr);
  if plotIndex ~= 3
    text(0.05, 0.90, textStr, 'Units', 'normalized');
  else
    text(0.05, 0.10, textStr, 'Units', 'normalized');
  end
  
  % plot a histogram for grand kerne, mark the CI
  limits = get(gca, 'YLim');
  limit = max(abs(limits));
  edges = linspace(-limit, limit, 12);
  subplot(2, 4, plotIndex + 4);
  h = histogram(grandKernel, edges);
  hold on;
  limits = get(gca, 'YLim');
  plot([-CI, -CI], limits, 'k:');
  plot([CI, CI], limits, 'k:');
  q = quantile(h.Data, [0.05, 0.95]);
  plot([q(1), q(1)], limits, 'r:');
  plot([q(2), q(2)], limits, 'r:');

end