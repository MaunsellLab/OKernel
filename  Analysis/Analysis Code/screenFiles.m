function screenFiles
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
  
	mode = 'normal';    % standard plots
% 	mode = 'control';     % control sessions with offset fiber
% 	mode = 'prePostControl';     % before and after control sessions
  
  dataDirName = '/Users/Shared/Data/OKernel/';
	tableDataName = [dataDirName ' Analysis/Processed Files.mat'];
  minLimits = 0;
  oneDay = [];
  switch mode
    case {'normal', 'test'}
      rampLimits = [0, 500];
      decLimits = 0;  
    otherwise
      rampLimits = 0;
      decLimits = -1;
      oneDay = [];
  end
  
% All animals, step and ramp
	animals = {'All'};
  
% Performance of individual step animals
%   rampLimits = 0;
% 	animals = {'844', '866', '902', '903', '905', '1003', '1112', '1145', '1150', '1218', '1220', '1223', '1257'};

% Performance of individual ramp animals (not used in a figure)
%     rampLimits = 500;
%   	animals = {'902', '1112', '1145', '1150', '1181', '1218', '1220', '1257'};

% Example session
%   animals = {'902'};
%   oneDay = '2019-10-10';

% Set up to plot the selected sessions
  switch mode
  case 'test'
    RTMinMS = 200; RTMaxMS = 500; minTrials = 0; minDecrement = 0;
    peakCutoff = 0; rampMS = 0;
    tableDataName = [dataDirName ' Analysis/Test Files.mat'];            
    [U, ~] = getSubset(mode, dataDirName, tableDataName, animals, rampMS, minTrials, minDecrement, oneDay);
    if size(U, 1) == 0
        fprintf('No valid sessions found for rampMS %d and threshold factor %.2f\n', rampMS, peakCutoff);
    else
        headerText = 'Test mode';
        doOneFigure(U, dataDirName, RTMinMS, RTMaxMS, headerText);
    end
  case {'normal', 'control', 'prePostControl'}
    limits.minSessions = 0;                         % require at least 8 sessions for each animal
    for r = rampLimits
      for c = 0.0
        for a = 1:length(animals)
          for t = minLimits
            for d = decLimits
              limits.rampMS = r;
              limits.criterion = c;
              limits.animal = animals{a};
              limits.minTrials = t;
              limits.minDec = d;
              doOneCase(mode, dataDirName, tableDataName, oneDay, limits);
            end
          end
        end
      end
    end
  end
end

%%
function doOneCase(mode, dataDirName, tableDataName, oneDay, limits)

  [U, maxDate] = getSubset(mode, dataDirName, tableDataName, oneDay, limits);
  if size(U, 1) == 0
    fprintf('No valid sessions found for rampMS %d, animal %s minTrials %d and decrement %.2f\n', ...
      rampMS, animal, minTrials, minDecrement);
    return;
  end
  headerText = cell(1, 7);
  if limits.rampMS == 0
    headerText{1} = sprintf('Visual Stimulus Step');
  else
    headerText{1} = sprintf('Visual Stimulus Ramp %d ms', limits.rampMS);
  end
	[RTMinMS, RTMaxMS, missMinMS, stimStr] = getRTParams(limits.rampMS);
  headerText{2} = sprintf('Hit times \\geq%d and <%d ms', RTMinMS, RTMaxMS);
  headerText{3} = sprintf('Miss times \\geq%d ms', missMinMS);
  if strcmp(limits.animal, 'All')
    headerText{4} = sprintf('%d sessions from %d animals', size(U, 1), length(unique(U.animal)));
  else
    headerText{4} = sprintf('%d sessions from Animal %s', size(U, 1), limits.animal);
  end
  headerText{5} = sprintf('%d trial minimum for hits/misses', limits.minTrials);
  if limits.minDec == -1
    headerText{6} = 'No required decrease in hit rate with opto';
    decStr = ' No HitDec.pdf';
  else
    headerText{6} = sprintf('Opto hits >=%.0f%% below unstim hits', limits.minDec * 100.0);
    decStr = sprintf(' %2.0f%% Hit Dec.pdf', limits.minDec * 100);
  end
  doOneFigure(U, dataDirName, limits.rampMS, headerText);
%   dayStr = sprintf('%d', maxDate);
%   minStr = sprintf(' %d Trials', minTrials);
  switch mode
    case 'normal'
      modeStr = '';
    case 'control'
      modeStr = ' Control';
    case 'prePostControl'
      modeStr = ' PrePostControl';
  end
  saveas(gcf, ['/Users/Shared/Data/OKernel/ Analysis/Figures/', modeStr, stimStr, ' ', limits.animal, decStr]);
end

%%

function doOneFigure(U, dataDirName, rampMS, headerText)

  display(U);
  [plotStartMS, plotEndMS, plotRTStartMS] = plotLimits();
  
  h = figure(1);
  set(h, 'Units', 'inches', 'Position', [25, 1.25, 8.5, 11]);
  clf;
  % Compile and plot the kernels
  % hit kernel
	offset = 0.5;
  ylabel = 'Normalized Power';
  % hit kernel
  numHits = sum(U.hits);
  plotTitle = sprintf('Hit Kernel (n=%d)', numHits);
  hitCI = stimCI(numHits);
  hitKernel = sum(cell2mat(U.hitKernel) .* double(U.hits), 1) / numHits;
  doOneKernelPlot(4, hitKernel, 'stim', plotStartMS, plotEndMS, plotTitle, ylabel, offset + hitCI, offset - hitCI);

  % miss kernel
  numMisses = sum(U.misses);
  plotTitle = sprintf('Miss Kernel (n=%d)', numMisses);
  missCI = stimCI(numMisses);
  missKernel = sum(cell2mat(U.missKernel) .* double(U.misses), 1) / numMisses;
  doOneKernelPlot(5, missKernel, 'stim', plotStartMS, plotEndMS, plotTitle, '', offset + missCI, offset - hitCI);

  % total kernel trials weighted across all trials. We need to multiple the weighted sum by 2 because it is effectively
  % a mean of the hit and miss kernels, not a difference. By taking the mean, we lose the doubling that we should get
  % from the opposing effects.  This has been validated in simulations. 
  plotTitle = sprintf('Weight by Trial (n=%d)', numHits + numMisses);
  kernel = (sum((cell2mat(U.hitKernel) - 0.5) .* double(U.hits), 1) - ...
      sum((cell2mat(U.missKernel) - 0.5) .* double(U.misses), 1)) / sum(U.hits + U.misses) * 2;
  totalCI = stimCI(sum(U.hits + U.misses));
  kernel = doOneKernelPlot(6, kernel, 'stim', plotStartMS, plotEndMS, plotTitle, '', totalCI, -totalCI);
  sigmas = min(kernel) / -totalCI * 1.96;                   % use 95% CI to get SEMs of kernel peak
  
  % RT aligned kernel
  plotTitle = sprintf('RT Aligned (n=%d)', numHits);
  doOneKernelPlot(7, mean(cell2mat(U.RTKernel), 1), 'RT', plotRTStartMS, plotRTStartMS + plotEndMS - plotStartMS, ...
      plotTitle, ylabel, offset + hitCI, offset - hitCI);

  % Stim-RT aligned kernel
  plotTitle = sprintf('Stim-RT Aligned (n=%d)', numHits);
  doOneKernelPlot(10, mean(cell2mat(U.SRTKernel), 1), 'stimRT', plotRTStartMS, plotRTStartMS + plotEndMS - plotStartMS, ...
      plotTitle, ylabel, offset + hitCI, offset - hitCI);

  % FA kernel. There might be some sessions with no FA, so we must clear them out
  kernelLength = zeros(1, size(U, 1));
  for index = 1:size(U, 1)
      kernelLength(index) = size(U.FAKernel{index}, 2);
  end
  validIndices = find(kernelLength > 1);
  if ~isempty(validIndices)
      FAKernel = U.FAKernel(validIndices);
      FAs = sum(U.FAs(validIndices));
      plotTitle = sprintf('FA aligned (n=%d)', sum(FAs));
      faCI = stimCI(FAs);
      doOneKernelPlot(8, mean(cell2mat(FAKernel), 1), 'FA', plotStartMS, plotEndMS, plotTitle, '', ...
              offset + faCI, offset - faCI);
  end

  % Random Kernel
  if numHits > 0 && numMisses > 0
      plotTitle = sprintf('Random Kernel (n=%d)', numHits + numMisses);
      doOneKernelPlot(9, mean(cell2mat(U.randomKernel), 1), 'stim', plotStartMS, plotEndMS, plotTitle, '', ...
        totalCI, -totalCI);
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
  [RTMinMS, RTMaxMS, ~, ~] = getRTParams(rampMS);
  doRTHistogramPlot(correctRTs, wrongRTs, missRTs, RTMinMS, RTMaxMS, minRespTimeMS, maxRespTimeMS);
  doRTPDFPlot(correctRTs, wrongRTs, missRTs, RTMinMS, RTMaxMS, minRespTimeMS, maxRespTimeMS)

  % Coordinate the scaling across plots
  sameYAxisScaling(4, 3, [4, 5, 7, 8, 10]);
  sameYAxisScaling(4, 3, [6, 9]);
%   sameYAxisScaling(4, 3, [6, 9], [-0.125, 0.075]);      % for plotting individual ramps
  
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
  axisHandle = subplot(4, 3, 1);						% default axes are 0 to 1
  set(axisHandle, 'visible', 'off');
  set(axisHandle, 'outerPosition', [0.02 0.75, 0.25, 0.2]);
  text(0.00, 1.25, 'OKernel', 'FontWeight', 'bold', 'FontSize', 16);
  headerText{end + 1} = sprintf('Hit rate no stim %.3f (SE %.3f)\n                  stim %.3f (SE %.3f)', ...
    noStimHitRate, noStimHitRateSE, stimHitRate, stimHitRateSE);
  headerText{end + 1} = sprintf('FA rate no stim %.3f (SE %.3f)\n                  stim %.3f (SE %.3f)', ...
    noStimFARate, noStimFARateSE, stimFARate, stimFARateSE);
  headerText{end + 1} = sprintf('Kernel peak at %.1f sigma', sigmas);
  text(0.00, 1.10, headerText, 'VerticalAlignment', 'top');

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
