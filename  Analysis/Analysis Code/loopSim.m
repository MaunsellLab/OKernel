function loopSim()

  hitRates = 0.20:0.05:1.00;
  nRates = length(hitRates);

  hitMissErr = zeros(1, nRates);
  oneKernelErr = zeros(1, nRates);
  weightedErr = zeros(1, nRates);
  pHit = zeros(1, nRates);

  for i = 1:nRates
    [pHit(i), hitMissErr(i), oneKernelErr(i), weightedErr(i)] = whiteNoiseOptoSim2(25000, false, hitRates(i));
  end
  h = figure(2);
  set(h, 'Units', 'inches', 'Position', [25, 1, 8.5, 6]);
  clf;
  plot(pHit, hitMissErr, pHit, oneKernelErr, pHit, weightedErr);
  legend('Hit/Miss', 'Single Kernel x2', 'Weighted Hit/Miss', 'location', 'south');
  ylabel('explained variance');
  ylim([0, 100]);
  xlim([0, 1.0]);
  xlabel('hit probability on opto trials');

end
