function [hitMissErr, singleErr, weightedErr] = whiteNoiseOptoSim(nTrials, silent)
% [hitMissErr, singleErr, weightedErr] = whiteNoiseOptoSim(nTrials [, silent])


%% Parameters

baseHitRate = 0.9;  % Hit rate on non-opto trials
pOpto = 0.5;        % Probability of a trial being an opto trial

yLims = [0.3 0.7];  % y-axis limits for plotting


%% Optional argument

if nargin < 2
  silent = 0;
end


%% Construct an opto influence timecourse

optoEffect = zeros(1, 20);
optoEffect(3:8) = -[1 3 3 2 1.5 0.5] / 20;


%% Construct random opto patterns

oStims = round(rand(nTrials, length(optoEffect)));


%% Decide which trials will get opto stim

oTrials = (rand(nTrials, 1) < pOpto);
nOpto = sum(oTrials);

%% Determine the outcome of each trial

hits = rand(nTrials, 1) < baseHitRate + oTrials .* sum(optoEffect .* oStims, 2);


%% Basic verification the sim worked

if ~silent
  fprintf('p(Hit) baseline: %0.1f%%\n', 100 * mean(hits(~oTrials)));
  fprintf('p(Hit) opto: %0.1f%%\n', 100 * mean(hits(oTrials)));
end


%% Overall kernel computed using hit kernel - miss kernel

hitKernel = mean(oStims(hits & oTrials, :));
missKernel = mean(oStims(~hits & oTrials, :));
hitMissKernel = hitKernel - missKernel;

hitMissErr = sum(abs(hitMissKernel - optoEffect));

if ~silent
  figure;
  
  subplot(1, 3, 1);
  hold on;
  plot(hitKernel, 'LineWidth', 2);
  plot(0.5 + optoEffect, 'r-', 'LineWidth', 1);
  ylim(yLims);
  title('Hit kernel')
  
  subplot(1, 3, 2);
  hold on;
  plot(missKernel, 'LineWidth', 2);
  plot(0.5 - optoEffect, 'r-', 'LineWidth', 1);
  ylim(yLims);
  title('Miss kernel');
  
  subplot(1, 3, 3);
  hold on;
  plot(hitMissKernel, 'LineWidth', 2);
  plot(optoEffect, 'r-', 'LineWidth', 1);
  ylim(yLims - 0.5);
  title('Hit - miss kernel');
  
  fprintf('Hit - Miss kernel error: %0.2f\n', hitMissErr);
end



%% Overall kernel computed using hit and -miss trials in one kernel

oneKernel = (sum(oStims(hits & oTrials, :)) + sum(1 - oStims(~hits & oTrials, :))) / nOpto - 0.5;

singleErr = sum(abs(oneKernel * 2 - optoEffect));

if ~silent
  figure;
  subplot(1, 2, 1);
  hold on;
  plot(oneKernel, 'LineWidth', 2);
  plot(optoEffect, 'r-', 'LineWidth', 1);
  ylim(yLims - 0.5);
  title('Single kernel');
  
  subplot(1, 2, 2);
  hold on;
  plot(2 * oneKernel, 'LineWidth', 2);
  plot(optoEffect, 'r-', 'LineWidth', 1);
  ylim(yLims - 0.5);
  title('Single kernel * 2');
  
  fprintf('Single kernel error: %0.2f\n', sum(abs(oneKernel - optoEffect)));
  fprintf('Single kernel * 2 error: %0.2f\n', singleErr);
end



%% Overall kernel computed using hit kernel - miss kernel with weighted averaging

nOpto = sum(oTrials);
nOptoHits = sum(hits & oTrials);
nOptoMisses = sum(~hits & oTrials);

hitWeight = nOpto / nOptoMisses;
hitKernel =  0.5 + hitWeight * (mean(oStims(hits & oTrials, :)) - 0.5);

missWeight = nOpto / nOptoHits;
missKernel = 0.5 + missWeight * (mean(oStims(~hits & oTrials, :)) - 0.5);

hitMissKernelWeighted = ((hitKernel - 0.5) * nOptoHits - (missKernel - 0.5) * nOptoMisses) / nOpto;

weightedErr = sum(abs(hitMissKernelWeighted - optoEffect));

if ~silent
  figure;
  
  subplot(1, 3, 1);
  hold on;
  plot(hitKernel, 'LineWidth', 2);
  plot(0.5 + optoEffect, 'r-', 'LineWidth', 1);
  ylim(yLims);
  title('Re-weighted Hit kernel')
  
  subplot(1, 3, 2);
  hold on;
  plot(missKernel, 'LineWidth', 2);
  plot(0.5 - optoEffect, 'r-', 'LineWidth', 1);
  ylim(yLims);
  title('Re-weighted Miss kernel');
  
  subplot(1, 3, 3);
  hold on;
  plot(hitMissKernelWeighted, 'LineWidth', 2);
  plot(optoEffect, 'r-', 'LineWidth', 1);
  ylim(yLims - 0.5);
  title('Weighted Total kernel');
  
  fprintf('Weighted kernel error: %0.2f\n', weightedErr);
end
