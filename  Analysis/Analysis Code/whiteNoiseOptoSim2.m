function [pOptoHit, hitMissErr, oneKernel2Err, weightedErr] = whiteNoiseOptoSim2(nTrials, silent, baseHitRate)
% [hitMissErr, weightedErr] = whiteNoiseOptoSim2(nTrials[, silent])

  %% Parameters

  pOpto = 0.5;        % Probability of a trial being an opto trial
  yLims = [0.1 0.9];  % y-axis limits for plotting

  %% Optional argument

  if nargin < 2
    silent = false;
  end
  if nargin < 3
    baseHitRate = 0.9;  % Hit rate on non-opto trials
  end

  %% Construct an opto influence timecourse

  optoEffect = zeros(1, 20);
  optoEffect(3:10) = -[3 3 3 2.5 2.0 1.5 1.0 0.5] / 50;

  %% Construct random opto patterns and assign hits

  oStims = round(rand(nTrials, length(optoEffect)));
  oTrials = (rand(nTrials, 1) < pOpto);
  nOpto = sum(oTrials);
  hits = rand(nTrials, 1) < baseHitRate + oTrials .* sum(optoEffect .* oStims, 2);
  nOptoHits = sum(hits & oTrials);
  nOptoMisses = sum(~hits & oTrials);
  pOptoHit = mean(hits(oTrials));

  %% Display header data

  if ~silent
    h = figure(1);
    set(h, 'Units', 'inches', 'Position', [25, 9, 8.5, 11]);
    clf;
    set(gca, 'visible', 'off');
    subplot(5, 3, 1);
    headerStr{1} = sprintf('%d trials, %d opto', nTrials, nOpto);
    headerStr{2} = sprintf('p(Hit) baseline: %0.1f%%', 100 * mean(hits(~oTrials)));
    headerStr{3} = sprintf('p(Hit) opto: %0.1f%%', 100 * pOptoHit);
    text(0.0, 0.95, headerStr, 'units', 'normalized', 'verticalAlignment', 'top', 'fontSize', 12);
    set(gca, 'visible', 'off');

  end

  %% Overall kernel computed using hit kernel - miss kernel

  hitKernel = mean(oStims(hits & oTrials, :));
  missKernel = mean(oStims(~hits & oTrials, :));
  hitMissKernel = hitKernel - missKernel;
  hitMissErr = explainedVarPC(optoEffect, hitMissKernel);

  if ~silent
    doOnePlot(1, hitKernel, yLims, 0.5 + optoEffect, 'Hit Kernel', nOptoHits);
    doOnePlot(2, missKernel, yLims, 0.5 - optoEffect, 'Miss Kernel', nOptoMisses);
    doOnePlot(3, hitMissKernel, yLims - 0.5, optoEffect, 'Hit - Miss Kernel', nOpto, hitMissErr);
  end

  %% Overall kernel computed using hit and -miss trials in one kernel

  oneKernel = (sum(oStims(hits & oTrials, :)) + sum(1 - oStims(~hits & oTrials, :))) / nOpto - 0.5;
  oneKernelErr = explainedVarPC(optoEffect, oneKernel);
  oneKernel2Err = explainedVarPC(optoEffect, oneKernel * 2);
  if ~silent
    doOnePlot(4, oneKernel, yLims - 0.5, optoEffect, 'Single Kernel', nOpto, oneKernelErr);
    doOnePlot(5, 2 * oneKernel, yLims - 0.5, optoEffect, 'Single kernel * 2', nOpto, oneKernel2Err);
  end

  %% Overall kernel computed using hit kernel - miss kernel with weighted averaging

  hitWeight = nOpto / nOptoMisses;
  hitKernel =  0.5 + hitWeight * (mean(oStims(hits & oTrials, :)) - 0.5);
  hitWeightedErr = explainedVarPC(optoEffect, hitKernel - 0.5);
  missWeight = nOpto / nOptoHits;
  missKernel = 0.5 + missWeight * (mean(oStims(~hits & oTrials, :)) - 0.5);
  missWeightedErr = explainedVarPC(optoEffect, -(missKernel - 0.5));
  hitMissKernelWeighted = ((hitKernel - 0.5) * nOptoHits - (missKernel - 0.5) * nOptoMisses) / nOpto;
  weightedErr = explainedVarPC(optoEffect, hitMissKernelWeighted);

  if ~silent
    doOnePlot(7, hitKernel, yLims, 0.5 + optoEffect, 'Re-weighted Hit Kernel', nOptoHits, hitWeightedErr);
    doOnePlot(8, missKernel, yLims, 0.5 - optoEffect, 'Re-weighted Miss Kernel', nOptoMisses, missWeightedErr);
    doOnePlot(9, hitMissKernelWeighted, yLims - 0.5, optoEffect, 'Weighted Total Kernel', nOpto, weightedErr);
  end
end

function doOnePlot(plotIndex, data, yLims, optoEffect, plotTitle, nTrials, error)

  subplot(5, 3, plotIndex + 3);
  hold on;
  plot(data, 'LineWidth', 2);
  plot(optoEffect, 'r-', 'LineWidth', 1);
  ylim(yLims);
  title(plotTitle);
	if nargin > 5
    if nargin > 6
      text(0.50, 0.95, sprintf('%d trials\n%0.0f%% variance', nTrials, error), ...
            'units', 'normalized', 'verticalAlignment', 'top');
    else
      text(0.55, 0.95, sprintf('%d trials', nTrials), 'units', 'normalized', 'verticalAlignment', 'top');
    end
  end
end

function eVarPC = explainedVarPC(y, yFit)

eVarPC = max(0.0, 100.0 * (1.0 - (mean((y - yFit).^2) / var(y))));

end

