function profile = doOneKernelPlot(subIndex, profile, type, startTimeMS, endTimeMS, plotTitle, yLabel, posCI, negCI)

  subplot(4, 3, subIndex);
  cla reset;
  bins = size(profile, 2);
  if bins < 25
      return;
  end
  if posCI ~= 0 && negCI ~= 0
      h = fill([0, bins, bins, 0], [posCI, posCI, negCI, negCI], [0.8, 0.8, 0.8]);
      set(h, 'linestyle', ':', 'facealpha', 0.25)
  end
  hold on;

  % low-pass filtering was worked out using testFiltering.m
  sampleFreqHz = 1000;
  filterLP = designfilt('lowpassfir', 'PassbandFrequency', 90 / sampleFreqHz, ...
    'StopbandFrequency', 2 * 90 / sampleFreqHz, 'PassbandRipple', 1, 'StopbandAttenuation', 60, ...
    'DesignMethod','equiripple');

  meanProfile = mean(profile);
  profile = filtfilt(filterLP, profile - meanProfile) + meanProfile;
  plot(profile, 'b');

  ax = gca;
  ax.XGrid = 'on';
  xlim(ax, [0, bins]);
  [plotStartMS, plotEndMS, plotRTStartMS] = plotLimits(true);
  [stimToRTMS, postRTMS] = stimRTLimits();
  preStimMS = -((endTimeMS - startTimeMS) - stimToRTMS - postRTMS);     % start relative to stim on
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
      set(gca,'XTick', [0, -preStimMS, -preStimMS + stimToRTMS, bins]);
      set(gca, 'XTickLabel', {sprintf('%d', preStimMS), 'Stim', 'RT', sprintf('+%d', endTimeMS)});
      xlabel('Fixed Stim-RT Interval');
    otherwise
      fprintf('doOneKernelPlot: unrecognized plot type: %s\n', type);
  end
  ylabel(yLabel);
  title(plotTitle);
  hold off;
end
