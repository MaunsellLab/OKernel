function plotSubset
  % Plot a page full of kernels for subset of the data
 
  [T, dataDirName, limits] = getSessionTable('All');

  if length(animals) > 1
    limits.fig = figure(2);
    set(limits.fig, 'Units', 'inches', 'Position', [25, 11.25, 8.5, 11]);
    clf;
    for s = 1:length(animals)
      limits.ax(s) = subplot(4, 3, s);
      limits.ax(s).Visible = 'Off';
    end
  end

  for r = rampLimits
    for a = 1:length(animals)
      limits.rampMS = r;
      limits.animal = animals{a};
      limits.aniNum = a;
      dataName = sprintf('Ramp %d%s', limits.rampMS);
      [U, ~] = getSubset(mode, dataDirName, [dataDirName, ' Analysis/Mat Files/', dataName, '.mat'], limits);
      if size(U, 1) == 0
        return;
      end
      bootstraps = getCaseBootstraps(T, dataDirName, dataName, limits);
      doOneFigure(U, dataDirName, dataName, limits, bootstraps);
    end
  end

  if length(animals) > 1
    figure(2);
    sameYAxisScaling(4, 3, 1:length(animals));
    saveas(gcf, [dataDirName, ' Analysis/Figures/', sprintf('Ramp %d%s Individuals.pdf', limits.rampMS, modeStr)]);
  end
end

%%
function doOneFigure(U, dataDirName, dataName, limits, bootstraps)
  % Compile and plot the kernels
  display(U);
  [plotStartMS, plotEndMS, plotRTStartMS] = plotLimits();
  
  h = figure(1);
  set(h, 'Units', 'inches', 'Position', [25, 1.25, 8.5, 11]);
  clf;
  ylabel = 'Normalized Power';
  limits.yAxis = 0.5;
  
  % hit kernel
  numHits = size(bootstraps.hitProfiles, 1);
  plotTitle = sprintf('Hit Kernel (n=%d)', numHits);
  subplot(4, 3, 4);
  CIs = doOneBootPlot(bootstraps.hitProfiles / 2 + 0.5, limits, 'stim', plotStartMS, plotEndMS, plotTitle, ylabel);
%   save(strcat(dataDirName, ' Analysis/Mat Files/', dataName, ' ', limits.animal, ' Hit Kernel'), 'CIs');
  
  % miss kernel
  numMisses = size(bootstraps.missProfiles, 1);
  plotTitle = sprintf('Miss Kernel (n=%d)', numMisses);
  subplot(4, 3, 5);
  doOneBootPlot(bootstraps.missProfiles / 2 + 0.5, limits, 'stim', plotStartMS, plotEndMS, plotTitle, '');

  % total kernel trials weighted across all trials. We need to multiple the weighted sum by 2 because it is effectively
  % a mean of the hit and miss kernels, not a difference. By taking the mean, we lose the doubling that we should get
  % from the opposing effects.  This has been validated in simulations. 
  plotTitle = sprintf('Weight by Trial (n=%d)', numHits + numMisses);
  limits.yAxis = 0.0;
  hitMissBoot = [bootstraps.hitProfiles; -bootstraps.missProfiles];
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
  doOneBootPlot(bootstraps.RTProfiles / 2 + 0.5, limits, 'RT', plotRTStartMS, ...
      plotRTStartMS + plotEndMS - plotStartMS, plotTitle, ylabel);
    
  % Stim-RT aligned kernel
  plotTitle = sprintf('Stim-RT Aligned (n=%d)', numHits);
  subplot(4, 3, 10);
  doOneBootPlot(bootstraps.stimRTProfiles / 2 + 0.5, limits, 'stimRT', plotRTStartMS, ...
      plotRTStartMS + plotEndMS - plotStartMS, plotTitle, ylabel);

  % FA kernel. There might be some sessions with no FA, so we must clear them out
	plotTitle = sprintf('FA Aligned (n=%d)', size(bootstraps.FAProfiles, 1));
	subplot(4, 3, 8);
  doOneBootPlot(bootstraps.FAProfiles / 2 + 0.5, limits, 'FA', plotStartMS, plotEndMS, plotTitle, '');

  % Random Kernel.  We randomly offset the hit/miss bootstraps to produce the random kernel
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
  minRespTimeMS = 10000; maxRespTimeMS = 0;
  rows = size(U, 1);
  for r = 1:rows
      clear file;
      load(strcat(dataDirName, U.animal(r), '/MatFiles/', U.date(r)), 'file');
      minRespTimeMS = min(minRespTimeMS, file.tooFastMS);         % set min/max response times
      maxRespTimeMS = max(maxRespTimeMS, file.rewardedLimitMS);
  end
  correctRTs = cat(2, U.correctRTs{:});
  earlyRTs = cat(2, U.earlyRTs{:});
  failRTs = cat(2, U.failRTs{:});
  doRTHistogramPlot(correctRTs, earlyRTs, failRTs, minRespTimeMS, maxRespTimeMS);
  doRTPDFPlot(correctRTs, earlyRTs, failRTs, minRespTimeMS, maxRespTimeMS)

  % Coordinate the scaling across plots
  sameYAxisScaling(4, 3, [4, 5, 7, 8, 10]);
  sameYAxisScaling(4, 3, [6, 9]);

  doHeader(U, limits);
  if ~isempty(limits.oneDay)
    saveas(gcf, strcat(dataDirName, ' Analysis/Figures/', dataName, ' ', limits.animal, ' ', limits.oneDay, '.pdf'));
  else
    saveas(gcf, strcat(dataDirName, ' Analysis/Figures/', dataName, ' ', limits.animal, '.pdf'));
  end
end