function testFiltering

% This function was used to look at the noise in average white noise
% kernels. It takes FFTs on random 25 ms kernels before and after
% filtering.  The upshot is that the kernels have a fair amount of noise
% above the pulse frequeuncy (25 ms or 40 Hz).  I used these plots to
% design a filter for excluding that noise.  It is best to use filtfilt()
% rather than filter(), because filter() messes up the first few entries in
% the kernel

  figure(1);
  clf;
  reps = [1, 100, 1000, 10000, 100000];
  numReps = length(reps);
  for r = 1:numReps
    doOnePlot(r, reps(r));
  end
  
  % do it again with a delta function to check whether that introduces higher frequencies
  figure(2);
  clf;
  reps = [1, 100, 1000, 10000, 100000];
  numReps = length(reps);
  for r = 1:numReps
    doOneDeltaPlot(r, reps(r));
  end

end

function doOneDeltaPlot(row, reps)

  sampleFreqHz = 1000;
  samplePeriodS = 1.0 / sampleFreqHz;
  numSamples = 5000;
%   times = (0:numSamples - 1) * samplePeriodS;
  pulseWidthMS = 25;
  pulseSamples = pulseWidthMS / (samplePeriodS * 1000);
  filterLP = designfilt('lowpassfir', 'PassbandFrequency', 60 / sampleFreqHz, ...
    'StopbandFrequency', 2 * 60 / sampleFreqHz, ...
    'PassbandRipple', 1, 'StopbandAttenuation', 60, 'DesignMethod','equiripple');
  
	subplot(5, 4, (row - 1) * 4 + 1);
  tSum = zeros(1, numSamples);
  for rep = 1:reps
    t = repelem(randi([0, 1], 1, numSamples / pulseSamples + 1), pulseSamples); 
    r = randi([0, pulseSamples - 1]);
    t = t(1 + r:numSamples + r);
    if t(numSamples / 2) == 1
      tSum = tSum + t;
    else
      tSum = tSum + (1 - t);
    end
  end
  kernel = tSum / reps;
  plot(kernel(numSamples / 2 - 500 + 1:numSamples /2  + 500));
  ylim([0.35, 0.65]);
  title(sprintf('%d trials', reps));
  
  n = 2^nextpow2(numSamples);
  Y = fft(kernel - mean(kernel), n);
  P2 = abs(Y / numSamples);
  P1 = P2(1:n / 2 + 1);
  P1(2:end-1) = 2 * P1(2:end-1);
  subplot(5, 4, (row - 1) * 4 + 2);
  plot(0:sampleFreqHz / n:(sampleFreqHz / 2 - sampleFreqHz/n), smooth(P1(1:n/2), 50))
  title('frequency domain');
  xlim([0, 100]);
  
  subplot(5, 4, (row - 1) * 4 + 3);
	kernel = filtfilt(filterLP, kernel);
  plot(kernel(1:1000));
  ylim([0.35, 0.65]);
  title(sprintf('%d trials filtered', reps));

	Y = fft(kernel - mean(kernel), n);
  P2 = abs(Y / numSamples);
  P1 = P2(1:n / 2 + 1);
  P1(2:end-1) = 2 * P1(2:end-1);
  subplot(5, 4, (row - 1) * 4 + 4);
  plot(0:sampleFreqHz / n:(sampleFreqHz / 2 - sampleFreqHz/n), smooth(P1(1:n/2), 50))
  title('frequency domain filtered kernel');
  xlim([0, 100]);
  

end

function doOnePlot(row, reps)

  sampleFreqHz = 1000;
  samplePeriodS = 1.0 / sampleFreqHz;
  numSamples = 5000;
%   times = (0:numSamples - 1) * samplePeriodS;
  pulseWidthMS = 25;
  pulseSamples = pulseWidthMS / (samplePeriodS * 1000);
  filterLP = designfilt('lowpassfir', 'PassbandFrequency', 60 / sampleFreqHz, ...
    'StopbandFrequency', 2 * 60 / sampleFreqHz, ...
    'PassbandRipple', 1, 'StopbandAttenuation', 60, 'DesignMethod','equiripple');
  
	subplot(5, 4, (row - 1) * 4 + 1);
  tSum = zeros(1, numSamples);
  for rep = 1:reps
    t = repelem(randi([0, 1], 1, numSamples / pulseSamples + 1), pulseSamples); 
    r = randi([0, pulseSamples - 1]);
    tSum = tSum + t(1 + r:numSamples + r);
  end
  kernel = tSum / reps;
  plot(kernel(1:1000));
  ylim([0.35, 0.65]);
  title(sprintf('%d trials', reps));
  
  n = 2^nextpow2(numSamples);
  Y = fft(kernel - mean(kernel), n);
  P2 = abs(Y / numSamples);
  P1 = P2(1:n / 2 + 1);
  P1(2:end-1) = 2 * P1(2:end-1);
  subplot(5, 4, (row - 1) * 4 + 2);
  plot(0:sampleFreqHz / n:(sampleFreqHz / 2 - sampleFreqHz/n), smooth(P1(1:n/2), 50))
  title('frequency domain');
  xlim([0, 100]);
  
  subplot(5, 4, (row - 1) * 4 + 3);
	kernel = filtfilt(filterLP, kernel);
  plot(kernel(1:1000));
  ylim([0.35, 0.65]);
  title(sprintf('%d trials filtered', reps));

	Y = fft(kernel - mean(kernel), n);
  P2 = abs(Y / numSamples);
  P1 = P2(1:n / 2 + 1);
  P1(2:end-1) = 2 * P1(2:end-1);
  subplot(5, 4, (row - 1) * 4 + 4);
  plot(0:sampleFreqHz / n:(sampleFreqHz / 2 - sampleFreqHz/n), smooth(P1(1:n/2), 50))
  title('frequency domain filtered kernel');
  xlim([0, 100]);
  

end
