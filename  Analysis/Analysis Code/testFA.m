function testFA

minMS = 500;          % minimum prestim time
maxMS = 3000;         % maximum prestim time
tooFastMS = 100;
pH = 0.5;             % probability of a hit
repsPerLambda = 25;  % number of repetitions
numTrials = 100;
lambdas = [0.005, 0.001, 0.0005, 0.0001];        % rate of 1/s for a time period of 1 ms

for l = 1:length(lambdas)
  lambda = lambdas(l);
  fitL = zeros(1, repsPerLambda);
  uFitL = zeros(1, repsPerLambda);
  numF = zeros(1, repsPerLambda);
  for r = 1:repsPerLambda
    [fitL(r), uFitL(r), numF(r)] = doOneFit(lambda, numTrials, minMS, maxMS, tooFastMS, pH);
  end
  fprintf('L = %.5f\n', lambda);
  fprintf(' Uncorrected L fit = %.5f SEM %.5f (average Fs = %.1f)\n', mean(uFitL), std(uFitL)/sqrt(repsPerLambda), mean(numF));
  fprintf(' Corrected   L fit = %.5f SEM %.5f (average Fs = %.1f)\n', mean(fitL), std(fitL)/sqrt(repsPerLambda), mean(numF));
end

% dprime(pH, numF / numTrials)

end

%%
%{
Put the Fs into n (5?) equally populated bins (to control for variance and emphasize the larger points. 
These bins could be corrected for the decline in emphasis, but they will be unevenly weighted.  Presumably this
wouldn't be too big a factor.  It might be possible to create a weighting vector with one entry per F and take
the average of this for each bin.

%}

function [fitL, uFitL, numF] = doOneFit(lambda, numTrials, minMS, maxMS, tooFastMS, pH)

numF = 0;
numH = 0;
pF = lambda * exp(-lambda); % Poisson probability for (1) F

delaysMS = [];
stimMSs = zeros(1, numTrials);
for n = 1:numTrials
  stimMS = minMS + (n - 1) * (maxMS - minMS) / (numTrials - 1);
  stimMSs(n) = stimMS;
  f = false;
  for t = 1:stimMS + tooFastMS
    if rand() < pF
      delaysMS = [delaysMS, t]; %#ok<*AGROW>
      numF = numF + 1;
      f = true;
      break;
    end
  end
  if ~f
    if rand() < pH
      numH = numH + 1;
    end
  end
end

dT = 250;
edges = 0:dT:maxMS;
trialWeights = zeros(1, length(edges) - 1);
b = 2;
while edges(b) < minMS                        % t < minMS always contributes to possible Fs
  trialWeights(b - 1) = 1.0;
  b = b + 1;
end
% give the bin straddling minMS a fraction of full trials, and a fraction of declining trials
trialWeights(b - 1) = (minMS - edges(b - 1)) / dT;
binMiddle = edges(b) - (edges(b) - minMS) / 2;
trialWeights(b - 1) = trialWeights(b - 1) +(1.0 - binMiddle * (1.0 - 1.0 / numTrials) / (maxMS - minMS)) * (edges(b) - minMS) / dT;
for b = b + 1:length(edges)
  binMiddle = edges(b) - dT / 2;
  trialWeights(b - 1) = 1.0 - (binMiddle - minMS) * (1.0 - 1.0 / numTrials) / (maxMS - minMS);
end
y = histcounts(delaysMS, edges);
x = edges(1:end - 1) + dT / 2.0;
lastY = find(y > 0, 1, 'last' );
y = y(1:lastY)';
x = x(1:lastY)';
f = fit(x, y, 'exp1');
uFitL = -f.b;
% plot(f, x, y);
% xlim([0, maxMS]);
y = y ./ trialWeights(1:lastY)';
f = fit(x, y, 'exp1');
fitL = -f.b;

figure(2);
clf;
rows = 2;
subplot(rows, 1, 1);
plot(f, x, y);
xlim([0, maxMS]);

% fprintf('L = %.5f, L fit = %.5f (nF = %d)\n', lambda, -f.b, length(delaysMS));
end