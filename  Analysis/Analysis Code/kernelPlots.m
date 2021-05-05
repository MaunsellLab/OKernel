function kernelPlots
  % (based on screenFiles, modified to compute bootstrap CIs on the kernels)
  % initialize variables
  
  % Control sessions (done at the end of data collection) are not included
  % in the primary analysis
  % I tried setting a minimum number of trials, but it didn't improve the
  % signal to noise
  % I tried various delta-performance settings.  Requiring that performance
  % is not better with PV cell stimulation seems pretty good.  More
  % stringent requirement quickly reduced the data set. More liberal
  % criteria didn't help the signal.
  % I tried selecting files with significant differences in variance before
  % and after stimOn, but that didn't help. 
  % I also tried shifting the RTMin/RTMax window for the ramp stimulus, but
  % that also didn't help
  
	mode = 'normal';           	% standard plots
% 	mode = 'control';        	% control sessions with offset fiber
% 	mode = 'prePostControl'; 	% before and after control sessions

	limits.numBoot = 10;        % number of boot strap runs

  dataDirName = '/Users/Shared/Data/OKernel/';
  limits.oneDay = [];
  switch mode
    case {'normal'}
      rampLimits = [0, 500];
      limits.minSessions = 0;
    otherwise
      rampLimits = 0;
  end
  
% % All animals, step and ramp
% 	animals = {'All'};
%   
% Performance of individual step animals
%   rampLimits = 0;
% 	animals = {'902', '905', '1112', '1145', '1223'};

% Performance of individual ramp animals (not used in a figure)
    rampLimits = 500;
  	animals = {'902', '1112', '1150'};

% Example session
%   animals = {'902'};
%   limits.oneDay = '2019-10-10';

% Set up to plot the selected sessions
  limits.criterion = 0.0;                         % no criterion requirement on kernel significance
  limits.minTrials = 0;                           % no minimum number of trials
  switch mode
    case 'normal'
      modeStr = '';
      limits.minSessions = 5;                     % require at least n sessions for each animal
      limits.minDec = -1;                         % stim trials can't have better performance
      limits.minDPrime = -1;                      % minimum d' in the no stim condition
      limits.minAvgDeltaDPrime = 0.10;               % minimum effect of opto stimulation
      limits.maxMeanPowerMW = 0.25;               % maximum average power over sesions
    case 'control'
      modeStr = ' Control';
      limits.minSessions = 0;                    	% no minimum for control sessions
      limits.minDec = -1;                       	% no performance limit for control sessions
    case 'prePostControl'
      modeStr = ' PrePostControl';
      limits.minSessions = 0;                     % no minimum for control sessions
      limits.minDec = -1;                       	% no performance limit for control sessions
  end

  % If we're doing multiple animals, prepare a second page so that we can collect up all the total kernels and
  % display them together at the same scaling

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
      bootstraps = getCaseBootstraps(U, dataDirName, dataName, limits);
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