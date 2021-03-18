function testFA

minMS = 0;          % minimum prestim time
maxMS = 3000;         % maximum prestim time
tooFastMS = 0;
pH = 0.5;             % probability of a hit
repsPerLambda = 100;  % number of repetitions
numTrials = 1000;
lambdas = [0.005, 0.001, 0.0005];        % rate of 1/s for a time period of 1 ms

for l = 1:length(lambdas)
  lambda = lambdas(l);
  fitL = zeros(1, repsPerLambda);
%   unFitL = zeros(1, repsPerLambda);
  numF = zeros(1, repsPerLambda);
  for r = 1:repsPerLambda
    while numF(r) == 0
      [fitL(r), numF(r)] = doTruncFit(lambda, numTrials, minMS, maxMS, tooFastMS, mod(r,10) == 0);
    end
  end
  fprintf('L = %.5f\n', lambda);
%   fprintf(' Uncorrected L fit = %.5f SEM %.5f (average Fs = %.1f)\n', mean(unFitL), std(unFitL)/sqrt(repsPerLambda), mean(numF));
  fprintf(' Pure   L fit = %.5f SEM %.5f (average Fs = %.1f)\n', mean(fitL), std(fitL)/sqrt(repsPerLambda), mean(numF));
end

% dprime(pH, numF / numTrials)

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
  adjust = (minMS - binEdgesMS(b)) / diff(binEdgesMS(b:b + 1));           % portion before minMS
  fullFrac = 1.0 - mean([binEdgesMS(b + 1), minMS]) / (maxMS - minMS);   	% frac full trials in changing part of bin
  adjust = adjust + fullFrac * (binEdgesMS(b + 1) - minMS) / diff(binEdgesMS(b:b + 1));
  y(b) = y(b) / adjust;
  b = b + 1;
  
  % The remaining bins have a linearly changing number of short trials across their full width.
  for b = b:length(binEdgesMS) - 1
    fullFrac = 1.0 - mean(binEdgesMS(b:b + 1)) / (maxMS - minMS);         % average full trials across entire bin
      y(b) = y(b) / fullFrac;
  end
end

%%
function  delaysMS = padTrials(delaysMS, pF, minMS, maxMS)
%
% Adjust for the trials that ended early
%
  numTrials = sum(delaysMS >= minMS);
  for n = 1:numTrials
    startMS = minMS + (n - 1) * (maxMS - minMS) / (numTrials - 1);
    for t = startMS:maxMS
      if rand() < pF
        delaysMS = [delaysMS, t]; %#ok<*AGROW>
        break;
      end
    end
  end
end

%%
function [fitL, numF] = doTruncFit(lambda, numTrials, minMS, maxMS, tooFastMS, doPlots)
%
% Fits where the trials stop early, uniformly between minMS and maxMS
%
  % numH = 0;
  pF = lambda * exp(-lambda); % Poisson probability for (1) F

  delaysMS = [];
  stimMSs = zeros(1, numTrials);
  for n = 1:numTrials
    stimMS = minMS + (n - 1) * (maxMS - minMS) / (numTrials - 1);
    stimMSs(n) = stimMS;
    for t = 1:stimMSs(n)
      if rand() < pF
        delaysMS = [delaysMS, t]; %#ok<*AGROW>
        break;
      end
    end
  end
  
  dT = 100;
  numF = length(delaysMS);
  binEdgesMS = 0:dT:maxMS;
  x = binEdgesMS(1:end - 1) + dT / 2.0;
  y = histcounts(delaysMS, binEdgesMS);
  
  % adjust for the abbreviated trials
  y = adjustForShortTrials(y, binEdgesMS, minMS, maxMS);
  
  % fit() doesn't work well if there are a lot of y == 0 entries at the end of the sequence.  These are undersampled,
  % and they can't be corrected for the number of short trials. 
%   validBins = y > 0;
%   if sum(validBins) < 3
%     numF = 0;
%     fitL = 0;
%     return;
%   end
%   y = y(1:lastY)';
%   x = x(1:lastY)';
%   y = y(validBins)';
%   x = x(validBins)';
  lastY = find(y > 0, 1, 'last' );
  y = y(1:lastY)';
  x = x(1:lastY)';
  fitted = fit(x, y, 'exp1', 'startPoint', [0, 0]);
  fitL = -fitted.b;

  if doPlots
    figure(2);
    clf;
    rows = 2;
    subplot(rows, 1, 1);
    plot(fitted, x, y);
  end
end

%%
function [unfitL, fitL, numF] = doFit(lambda, numTrials, minMS, maxMS, tooFastMS, doPlots)

numF = 0;
% numH = 0;
numEdges = 11;
pF = lambda * exp(-lambda); % Poisson probability for (1) F

delaysMS = [];
stimMSs = zeros(1, numTrials);
for n = 1:numTrials
  stimMS = minMS + (n - 1) * (maxMS - minMS) / (numTrials - 1);
  stimMSs(n) = stimMS;
  stimMSs(n) = maxMS;
%   f = false;
  for t = 1:stimMS + tooFastMS
    if rand() < pF
      delaysMS = [delaysMS, t]; %#ok<*AGROW>
      numF = numF + 1;
%       f = true;
      break;
    end
  end
%   if ~f
%     if rand() < pH
%       numH = numH + 1;
%     end
%   end
end

binEdgesMS = prctile(delaysMS, 0:100 / (numEdges - 1):100);           % edges surrounding bins
binWidthsMS = diff(binEdgesMS);                                       % width of each bin
binMiddlesMS = binEdgesMS(1:numEdges - 1) + binWidthsMS / 2.0;        % middle of each bin
y = histcounts(delaysMS, binEdgesMS) ./ binWidthsMS * binWidthsMS(1); % counts adjusted for bin widths

weights = zeros(1, numEdges - 1);
for b = 1:numEdges - 1
  if binMiddlesMS(b) < minMS
    weights(b) = 1.0;
  else
    weights(b) = 1.0 - (binMiddlesMS(b) - minMS) * (1.0 - 1.0 / numTrials) / (maxMS - minMS);
  end
end
yCorrected = y ./ weights;

% dT = 250;
% edges = 0:dT:maxMS;
% trialWeights = zeros(1, length(edges) - 1);
% b = 2;
% while edges(b) < minMS                        % t < minMS always contributes to possible Fs
%   trialWeights(b - 1) = 1.0;
%   b = b + 1;
% end
% give the bin straddling minMS a fraction of full trials, and a fraction of declining trials
% trialWeights(b - 1) = (minMS - edges(b - 1)) / dT;
% binMiddle = edges(b) - (edges(b) - minMS) / 2;
% trialWeights(b - 1) = trialWeights(b - 1) +(1.0 - binMiddle * (1.0 - 1.0 / numTrials) / (maxMS - minMS)) * (edges(b) - minMS) / dT;
% for b = b + 1:length(edges)
%   binMiddle = edges(b) - dT / 2;
%   trialWeights(b - 1) = 1.0 - (binMiddle - minMS) * (1.0 - 1.0 / numTrials) / (maxMS - minMS);
% end
% y = histcounts(delaysMS, edges);
% x = edges(1:end - 1) + dT / 2.0;
% lastY = find(y > 0, 1, 'last' );
% y = y(1:lastY)';
% x = x(1:lastY)';
unfit = fit(binMiddlesMS', y', 'exp1');
unfitL = -unfit.b;
fitted = fit(binMiddlesMS', yCorrected', 'exp1');
fitL = -fitted.b;

% plot(f, x, y);
% xlim([0, maxMS]);
% y = y ./ trialWeights(1:lastY)';
% f = fit(x, y, 'exp1');
% fitL = -f.b;
if doPlots
  figure(2);
  clf;
  rows = 2;
  subplot(rows, 1, 1);
  plot(unfit, binMiddlesMS, y);
  xlim([0, maxMS]);
  hold on;
  % subplot(rows, 1, 2);
  plot(fitted, binMiddlesMS, yCorrected);
  % xlim([0, maxMS]);
end

end
%%
function [fitL, numF] = doPureFits(lambda, numTrials, minMS, maxMS, tooFastMS, doPlots)

  numF = 0;
  % numH = 0;
  pF = lambda * exp(-lambda); % Poisson probability for (1) F

  delaysMS = [];
  stimMSs = zeros(1, numTrials);
  for n = 1:numTrials
    stimMSs(n) = maxMS;
    for t = 1:stimMSs(n)
      if rand() < pF
        delaysMS = [delaysMS, t]; %#ok<*AGROW>
        numF = numF + 1;
        break;
      end
    end
  end

  dT = 100;
  binEdgesMS = 0:dT:maxMS;
  y = histcounts(delaysMS, binEdgesMS); % counts adjusted for bin widths
  x = binEdgesMS(1:end - 1) + dT / 2.0;
  lastY = find(y > 0, 1, 'last' );
  y = y(1:lastY)';
  x = x(1:lastY)';
  fitted = fit(x, y, 'exp1', 'startPoint', [0, 0]);
  fitL = -fitted.b;

  if doPlots
    figure(2);
    clf;
    rows = 2;
    subplot(rows, 1, 1);
    plot(fitted, x, y);
  end
end