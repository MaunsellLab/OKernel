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
% 	mode = 'control';         	% control sessions with offset fiber
% 	mode = 'prePostControl';      % before and after control sessions

	limits.numBoot = 10;     	% number of boot strap runs

  dataDirName = '/Users/Shared/Data/OKernel/';
  limits.oneDay = [];
  switch mode
    case {'normal', 'test'}
      rampLimits = [0];
      limits.minSessions = 0;
    otherwise
      rampLimits = 0;
  end
  
% All animals, step and ramp
	animals = {'All'};
  
% Performance of individual step animals
%   rampLimits = 0;
% 	animals = {'866', '902', '905', '1112', '1145', '1150', '1218', '1220', '1223', '1257'};

% Performance of individual ramp animals (not used in a figure)
%     rampLimits = 500;
%   	animals = {'902', '1112', '1145', '1150', '1181', '1218', '1220', '1257'};

% Example session
%   animals = {'902'};
%   limits.oneDay = '2019-10-10';

% Set up to plot the selected sessions
  switch mode
  case 'test'
    RTMinMS = 200; RTMaxMS = 500; minTrials = 0; minDecrement = 0;
    peakCutoff = 0; rampMS = 0;
    tableDataName = [dataDirName ' Analysis/Test Files.mat'];            
    [U, ~] = getSubset(mode, dataDirName, tableDataName, animals, rampMS, minTrials, minDecrement, limits.oneDay);
    if size(U, 1) == 0
        fprintf('No valid sessions found for rampMS %d and threshold factor %.2f\n', rampMS, peakCutoff);
    else
        headerText = 'Test mode';
        doOneFigure(U, dataDirName, RTMinMS, RTMaxMS, headerText);
    end
  case {'normal', 'control', 'prePostControl'}
   	limits.criterion = 0.0;                         % no criterion requirement on kernel significance
    limits.minTrials = 0;                           % no minimum number of trials
    switch mode
      case 'normal'
        modeStr = '';
        limits.minSessions = 8;                    	% require at least 8 sessions for each animal
        limits.minDec = 0;                              % stim trials can't have better performance
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
        doOneCase(mode, dataDirName, sprintf('Ramp %d%s', limits.rampMS, modeStr), limits);
      end
    end
    
    if length(animals) > 1
      figure(2);
      sameYAxisScaling(4, 3, 1:length(animals));
      saveas(gcf, [dataDirName, ' Analysis/Figures/', sprintf('Ramp %d%s Individuals.pdf', limits.rampMS, modeStr)]);
    end
  end
end

%%
function doOneCase(mode, dataDirName, dataName, limits)

  [U, ~] = getSubset(mode, dataDirName, [dataDirName, ' Analysis/Mat Files/', dataName, '.mat'], limits.oneDay, limits);
  if size(U, 1) == 0
    return;
  end
  bootstraps = getCaseBootstraps(U, dataDirName, dataName, limits);
  doOneFigure(U, dataDirName, dataName, limits, bootstraps);
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
  ax = subplot(4, 3, 4);
  doOneBootPlot(bootstraps.hitProfiles / 2 + 0.5, limits, 'stim', plotStartMS, plotEndMS, plotTitle, ylabel);
  
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
  wrongRTs = cat(2, U.wrongRTs{:});
  missRTs = cat(2, U.missRTs{:});
  [RTMinMS, RTMaxMS, ~, ~] = getRTParams(limits.rampMS);
  doRTHistogramPlot(correctRTs, wrongRTs, missRTs, RTMinMS, RTMaxMS, minRespTimeMS, maxRespTimeMS);
  doRTPDFPlot(correctRTs, wrongRTs, missRTs, RTMinMS, RTMaxMS, minRespTimeMS, maxRespTimeMS)

  % Coordinate the scaling across plots
  sameYAxisScaling(4, 3, [4, 5, 7, 8, 10]);
  sameYAxisScaling(4, 3, [6, 9]);
  
  % compute overall hit and FA rates
  
  numStim = sum(U.numStim);
  numStimHits = sum(U.hits);
  numStimMisses = sum(U.misses);
  stimHitRate = numStimHits / (numStimHits + numStimMisses);
  stimHitRateSE = sqrt(stimHitRate * (1 - stimHitRate) / (numStimHits + numStimMisses));
  stimFARate = sum(U.FAs) / numStim;
  stimFARateSE = sqrt(stimFARate * (1 - stimFARate) / numStim);

  numNoStim = sum(U.numNoStim);
  numNoStimHits = sum(U.noStimHits);
  numNoStimMisses = sum(U.noStimMisses);
  noStimHitRate = numNoStimHits / (numNoStimHits + numNoStimMisses);
  noStimHitRateSE = sqrt(noStimHitRate * (1 - noStimHitRate) / (numNoStimHits + numNoStimMisses));
  noStimFARate = sum(U.noStimFAs) / numNoStim;
  noStimFARateSE = sqrt(noStimFARate * (1 - noStimFARate) / numNoStim);
  
  % display header text
  
  headerText = cell(1, 7);
  if limits.rampMS == 0
    headerText{1} = sprintf('Visual Stimulus Step');
  else
    headerText{1} = sprintf('Visual Stimulus Ramp %d ms', limits.rampMS);
  end
	[RTMinMS, RTMaxMS, missMinMS, ~] = getRTParams(limits.rampMS);
  headerText{2} = sprintf('Hit times \\geq%d and <%d ms', RTMinMS, RTMaxMS);
  headerText{3} = sprintf('Miss times \\geq%d ms', missMinMS);
  if strcmp(limits.animal, 'All')
    headerText{4} = sprintf('%d sessions from %d animals', size(U, 1), length(unique(U.animal)));
  else
    headerText{4} = sprintf('%d sessions from Animal %s', size(U, 1), limits.animal);
  end
  headerText{5} = sprintf('%d bootstraps for CIs', limits.numBoot);
  if limits.minDec == -1
    headerText{6} = 'No required decrease in hit rate with opto';
  else
    headerText{6} = sprintf('Opto hits >=%.0f%% below unstim hits', limits.minDec * 100.0);
  end

  axisHandle = subplot(4, 3, 1);						% default axes are 0 to 1
  set(axisHandle, 'visible', 'off');
  set(axisHandle, 'outerPosition', [0.02 0.75, 0.25, 0.2]);
  text(0.00, 1.25, 'OKernel', 'FontWeight', 'bold', 'FontSize', 16);
  headerText{end + 1} = sprintf('Hit rate no stim %.3f (SE %.3f)\n                  stim %.3f (SE %.3f)', ...
    noStimHitRate, noStimHitRateSE, stimHitRate, stimHitRateSE);
  headerText{end + 1} = sprintf('FA rate no stim %.3f (SE %.3f)\n                  stim %.3f (SE %.3f)', ...
    noStimFARate, noStimFARateSE, stimFARate, stimFARateSE);
  text(0.00, 1.10, headerText, 'VerticalAlignment', 'top');
  if ~isempty(limits.oneDay)
    saveas(gcf, [dataDirName, ' Analysis/Figures/', dataName, ' ', limits.animal, ' ', limits.oneDay, '.pdf']);
  else
    saveas(gcf, [dataDirName, ' Analysis/Figures/', dataName, ' ', limits.animal, '.pdf']);
  end
end

%%
function [RTMinMS, RTMaxMS, missMinMS, stimStr] = getRTParams(rampMS)

  switch rampMS
    case 0
      RTMinMS = 200;
      RTMaxMS = 500;
      missMinMS = 750;
     	stimStr = 'Steps';
    case 500
      RTMinMS = 200;
      RTMaxMS = 500;
      missMinMS = 750;
      stimStr = 'Ramps';
    otherwise
      RTMinMS = 0;
      RTMaxMS = 0;
      missMinMS = 0;
      stimStr = 'Unknown';
  end
end
