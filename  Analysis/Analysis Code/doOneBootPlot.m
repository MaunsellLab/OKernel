function CIs = doOneBootPlot(bootstraps, limits, type, startTimeMS, endTimeMS, plotTitle, yLabel)

	sampleFreqHz = 1000;
  filterLP = designfilt('lowpassfir', 'PassbandFrequency', 90 / sampleFreqHz, ...
    'StopbandFrequency', 2 * 90 / sampleFreqHz, 'PassbandRipple', 1, 'StopbandAttenuation', 60, ...
    'DesignMethod','equiripple');

  bootKernel = bootstrp(limits.numBoot, @mean, bootstraps);
  PCs = prctile(bootKernel, [15.9, 50, 84.1]);              % +/- 1 SEM
  PCMeans = mean(PCs, 2);
  CIs = zeros(3, size(bootstraps, 2));
  for c = 1:3
    CIs(c, :) = filtfilt(filterLP, PCs(c, :) - PCMeans(c)) + PCMeans(c);
  end
  cla reset;
  bins = size(bootstraps, 2);
  if bins < 25
      return;
  end
  
	x = 1:size(CIs, 2);
  plot(x, CIs(2, :), 'b');
  hold on;
  x2 = [x, fliplr(x)];
  fillCI = [CIs(1, :), fliplr(CIs(3, :))];
  fill(x2, fillCI, 'b', 'lineStyle', '-', 'edgeColor', 'b', 'edgeAlpha', 0.5, 'faceAlpha', 0.10);
  ax = gca;
  xlim(ax, [0, bins]);
  ax.XGrid = 'on';
  [plotStartMS, plotEndMS, plotRTStartMS] = plotLimits();
	plot([0, bins], [limits.yAxis, limits.yAxis], 'k-');
  switch type
    case {'stim', 'Stim'}
      set(gca,'XTick', [0, -plotStartMS, -plotStartMS + 100, bins]);
      set(gca, 'XTickLabel', {sprintf('%d', plotStartMS), '0', '', sprintf('%d', plotEndMS)});
      xlabel('Time Relative to Stimulus');
    case {'rt', 'RT'}
      set(gca,'XTick', [0, -plotRTStartMS, bins]);
      set(gca, 'XTickLabel', {sprintf('%d', plotRTStartMS), '0', sprintf('%d', plotRTStartMS + plotEndMS - plotStartMS)});
      xlabel('Time Relative to RT');
    case {'early', 'Early'}
      set(gca,'XTick', [0, -plotRTStartMS, bins]);
      set(gca, 'XTickLabel', {sprintf('%d', plotRTStartMS), '0', sprintf('%d', plotRTStartMS + plotEndMS - plotStartMS)});
      xlabel('Time Relative to FA');
    case {'stimRT', 'StimRT'}
      [stimToRTMS, postRTMS] = stimRTLimits();
      preStimMS = -((endTimeMS - startTimeMS) - stimToRTMS - postRTMS);     % start relative to stim on
      set(gca,'XTick', [0, -preStimMS, -preStimMS + stimToRTMS, bins]);
      set(gca, 'XTickLabel', {sprintf('%d', preStimMS), 'Stim', 'RT', sprintf('+%d', endTimeMS)});
      xlabel('Fixed Stim-RT Interval');
    otherwise
      fprintf('doOneBootPlot: unrecognized plot type: %s\n', type);
  end
  ylabel(yLabel);
  title(plotTitle);
  hold off;
end
