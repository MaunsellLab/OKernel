function plotKernelPage(U, limits, stimProfiles)
  % Compile and plot the kernels
  display(U);
  [plotStartMS, plotEndMS, plotRTStartMS] = plotLimits();
  
  h = figure(1);
  set(h, 'Units', 'inches', 'Position', [25, 1.25, 8.5, 11]);
  clf;
  ylabel = 'Normalized Power';
  limits.yAxis = 0.5;
  limits.numBoot = 10;
  
  % hit kernel
  numHits = size(stimProfiles.hitProfiles, 1);
  plotTitle = sprintf('Hit Kernel (n=%d)', numHits);
  subplot(4, 3, 4);
  doOneBootPlot(stimProfiles.hitProfiles / 2 + 0.5, limits, 'stim', plotStartMS, plotEndMS, plotTitle, ylabel);
  
  % miss kernel
  numMisses = size(stimProfiles.missProfiles, 1);
  plotTitle = sprintf('Miss Kernel (n=%d)', numMisses);
  subplot(4, 3, 5);
  doOneBootPlot(stimProfiles.missProfiles / 2 + 0.5, limits, 'stim', plotStartMS, plotEndMS, plotTitle, '');

  % total kernel trials weighted across all trials. We need to multiple the weighted sum by 2 because it is effectively
  % a mean of the hit and miss kernels, not a difference. By taking the mean, we lose the doubling that we should get
  % from the opposing effects.  This has been validated in simulations. 
  plotTitle = sprintf('Weight by Trial (n=%d)', numHits + numMisses);
  limits.yAxis = 0.0;
  hitMissBoot = [stimProfiles.hitProfiles; -stimProfiles.missProfiles];
  subplot(4, 3, 6);
  doOneBootPlot(hitMissBoot, limits, 'stim', plotStartMS, plotEndMS, plotTitle, '');
  if ~strcmp(limits.animal, 'All')
    figure(2);
    h = subplot(4, 3, limits.aniNum);
    doOneBootPlot(hitMissBoot, limits, 'stim', plotStartMS, plotEndMS, plotTitle, '');
    h.Title.String = strrep(h.Title.String, 'Weight by Trial', ['Animal ', limits.animal]);
    figure(1);
  end
  
  % RT aligned kernel
  limits.yAxis = 0.5;
  plotTitle = sprintf('RT Aligned (n=%d)', numHits);
  subplot(4, 3, 7);
  doOneBootPlot(stimProfiles.RTProfiles / 2 + 0.5, limits, 'RT', plotRTStartMS, ...
      plotRTStartMS + plotEndMS - plotStartMS, plotTitle, ylabel);
    
  % Stim-RT aligned kernel
  plotTitle = sprintf('Stim-RT Aligned (n=%d)', numHits);
  subplot(4, 3, 10);
  doOneBootPlot(stimProfiles.stimRTProfiles / 2 + 0.5, limits, 'stimRT', plotRTStartMS, ...
      plotRTStartMS + plotEndMS - plotStartMS, plotTitle, ylabel);

  % Early kernel. There might be some sessions with no FA, so we must clear them out
	plotTitle = sprintf('Early Aligned (n=%d)', size(stimProfiles.earlyProfiles, 1));
	subplot(4, 3, 8);
  doOneBootPlot(stimProfiles.earlyProfiles / 2 + 0.5, limits, 'Early', plotStartMS, plotEndMS, plotTitle, '');

  % Random Kernel.  We randomly offset the hit/miss stimProfiles to produce the random kernel
  limits.yAxis = 0.0;
  if numHits > 0 && numMisses > 0
      arrayWidth = size(hitMissBoot, 2);
      for b = 1:size(hitMissBoot, 1)
        mult = randi([0 1], 1, 1) * 2 - 1;
        hitMissBoot(b, :) = circshift(hitMissBoot(b, :), randi([1, arrayWidth], 1, 1)) * mult;
      end
      plotTitle = sprintf('Random Kernel (n=%d)', numHits + numMisses);
      subplot(4, 3, 9);
      doOneBootPlot(hitMissBoot, limits, 'stim', plotStartMS, plotEndMS, plotTitle, '');
  end  
  
  %% Compile and plot the RT distributions
  minRespTimeMS = min(U.startRT(:));
  maxRespTimeMS = min(U.endRT(:));
  correctRTs = cat(2, U.correctRTs{:});
  earlyRTs = cat(2, U.earlyRTs{:});
  failRTs = cat(2, U.failRTs{:});
  doRTHistogramPlot(correctRTs, earlyRTs, failRTs, minRespTimeMS, maxRespTimeMS);
  doRTPDFPlot(correctRTs, earlyRTs, failRTs, minRespTimeMS, maxRespTimeMS)

  % Coordinate the scaling across plots
  sameYAxisScaling(4, 3, [4, 5, 7, 8, 10]);
  sameYAxisScaling(4, 3, [6, 9]);

  doHeader(U, limits);
end