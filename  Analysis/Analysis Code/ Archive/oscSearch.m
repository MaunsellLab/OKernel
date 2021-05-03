function oscSearch()
% Oscillation search

%{
The optogenetic kernels we are getting all seem to be dominated by an oscillation at 40 ms / 25 Hz. I have no idea
where this come from.  This simulation is to check whether it might be a natural consequence of the stimulus sequences
and filtering that we do.
%}

  % create kernel
  lenStim = 1000;
  pulseWidth = 25;                             % width of pulse (nominal ms)
  nTrials = 100000;                           % number of trials per session
  sampleFreqHz = 1000.0;

  bValues = randi([0, 1], nTrials, int32(lenStim / pulseWidth));    % random binary values (one per pulse)
  optoStim = repelem(bValues, 1, pulseWidth);                      	% expand to fill pulsewidths
  % give each stimulus a random phase by circularly shifting by a random portion of pulse width
  for t = 1:nTrials
    optoStim(t, :) = circshift(optoStim(t, :), randi([0, pulseWidth - 1]), 2);
  end
	kernel = mean(optoStim);
  
  k = zeros(1, lenStim);
  k(lenStim / 2 - pulseWidth:lenStim / 2) = 0:pulseWidth;
  k(lenStim / 2 + 1:lenStim / 2 + pulseWidth) = pulseWidth - 1:-1:0;
  k = k + rand(1, 1000) + pulseWidth / 5;
  kernel = k;
  
	filterLP = designfilt('lowpassfir', 'PassbandFrequency', 60 / sampleFreqHz, ...
    'StopbandFrequency', 2 * 60 / sampleFreqHz, 'PassbandRipple', 1, 'StopbandAttenuation', 60, ...
    'DesignMethod', 'equiripple');
  filtKernel = filtfilt(filterLP, kernel - 0.5);
	filterLP2 = designfilt('lowpassfir', 'PassbandFrequency', 12 / sampleFreqHz, ...
    'StopbandFrequency', 48 / sampleFreqHz, 'PassbandRipple', 1, 'StopbandAttenuation', 60, ...
    'DesignMethod', 'kaiserwin');
  filtKernel2 = filtfilt(filterLP2, kernel - 0.5);
  
  % plot some examples from the set of kernels
  % set up figure
  h = figure(1);
  set(h, 'Units', 'inches', 'Position', [25, 14.5, 10.24, 6.40]);
  clf;
  
  % plot the kernels and the CI
  CI = stimCI(nTrials);           % 95% CI for these kernels
  subplot(2, 3, 1);
  h = fill([0, size(kernel, 2), size(kernel, 2), 0], [CI, CI, -CI, -CI], [0.8, 0.8, 0.8]);
  set(h, 'linestyle', ':', 'facealpha', 0.25);
  hold on;
  plot(kernel - 0.5, 'b');
  xlim([0, size(kernel, 2)]);
  title('Raw Kernel');
  
  subplot(2, 3, 2);
  h = fill([0, size(kernel, 2), size(kernel, 2), 0], [CI, CI, -CI, -CI], [0.8, 0.8, 0.8]);
  set(h, 'linestyle', ':', 'facealpha', 0.25);
  hold on;
  plot(filtKernel, 'r');
  xlim([0, size(kernel, 2)]);
  title('LP Filter, 60 dB down between 60-120 Hz');
  
  subplot(2, 3, 3);
  h = fill([0, size(kernel, 2), size(kernel, 2), 0], [CI, CI, -CI, -CI], [0.8, 0.8, 0.8]);
  set(h, 'linestyle', ':', 'facealpha', 0.25);
  hold on;
  plot(filtKernel2, 'g');
  xlim([0, size(kernel, 2)]);
   title('LP Filter, 60 dB down between 12-48 Hz');
 
  subplot(2, 1, 2);
  n = 2^nextpow2(lenStim);
  Y = fft(kernel - mean(kernel), n);
  P2 = abs(Y / lenStim);
  P1 = P2(1:n / 2 + 1);
  P1(2:end-1) = 2 * P1(2:end-1);
  loglog(0:sampleFreqHz / n:(sampleFreqHz / 2 - sampleFreqHz/n), smooth(P1(1:n/2), 50), 'b');
  hold on;
  title('frequency domain');
  xlim([0, 100]);
  
  Y = fft(filtKernel - mean(filtKernel), n);
  P2 = abs(Y / lenStim);
  P1 = P2(1:n / 2 + 1);
  P1(2:end-1) = 2 * P1(2:end-1);
%   subplot(2, 2, 3);
  loglog(0:sampleFreqHz / n:(sampleFreqHz / 2 - sampleFreqHz/n), smooth(P1(1:n/2), 50), 'r');
  title('frequency domain');

  xlim([0, 100]);
  Y = fft(filtKernel2 - mean(filtKernel2), n);
  P2 = abs(Y / lenStim);
  P1 = P2(1:n / 2 + 1);
  P1(2:end-1) = 2 * P1(2:end-1);
%   subplot(2, 2, 3);
  loglog(0:sampleFreqHz / n:(sampleFreqHz / 2 - sampleFreqHz/n), smooth(P1(1:n/2), 50), 'g');
  title('frequency domain');
  xlim([1, 200]);

end