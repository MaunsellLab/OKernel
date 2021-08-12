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
  minLimits = 10;
  limits.oneDay = [];
  switch mode
    case {'normal'}
      rampLimits = [0];
      decLimits = 0.10;  
    otherwise
      rampLimits = 0;
      decLimits = -1;
  end
  
% All animals, step and ramp
	animals = {'All'};
  
% Performance of individual step animals
%   rampLimits = 0;
% 	animals = {'902', '905', '1150', '1145', '1112'};

% Performance of individual ramp animals (not used in a figure)
%     rampLimits = 500;
% 	animals = {'902', '1150', '1218', '1220'};

% Example session
%   animals = {'902'};
%   limits.oneDay = '2019-10-10';

% Set up to plot the selected sessions
  limits.minSessions = 10;                         % require at least n sessions for each animal
  limits.minDPrime = 0.5;
  for r = rampLimits
    limits.rampMS = r;
    for c = 0.0
    	limits.criterion = c;
      for a = 1:length(animals)
      	limits.animal = animals{a};
        for t = minLimits
          limits.minTrials = t;
          for d = decLimits
            limits.minDec = d;
            doCase(mode, dataDirName, tableDataName, limits);
          end
        end
      end
    end
  end
end

%%
function doCase(mode, dataDirName, tableDataName, limits)

  [U, ~] = getSubset(mode, dataDirName, tableDataName, limits);
  if size(U, 1) == 0
    fprintf('No valid sessions found for rampMS %d, animal %s minTrials %d and delta-d'' %.2f\n', ...
      limits.rampMS, limits.animal, limits.minTrials, limits.minDec);
    return;
  end
  doFigure(U, dataDirName, limits);
  switch mode
    case 'normal'
      modeStr = '';
    case 'control'
      modeStr = ' Control';
    case 'prePostControl'
      modeStr = ' PrePostControl';
  end
  if limits.rampMS == 0
    stimStr = 'Step';
  else
    stimStr = 'Ramp';
  end
  decStr = sprintf(' %.2f dP Dec', limits.minDec);
  saveas(gcf, ['/Users/Shared/Data/OKernel/ Analysis/Figures/', modeStr, stimStr, ' ', limits.animal, decStr '.pdf']);
end

%%

function doFigure(U, dataDirName, limits)

  display(U);
  [plotStartMS, plotEndMS, plotRTStartMS] = plotLimits();
  
  h = figure(1);
  set(h, 'Units', 'inches', 'Position', [25, 1.25, 8.5, 11]);
  clf;
  % Compile and plot the kernels
  % hit kernel
	offset = 0.5;
  ylabel = 'Normalized Power';
  
  % correct kernel
  numCorrect = sum(U.stimCorrects);
  plotTitle = sprintf('Hit Kernel (n=%d)', numCorrect);
  correctCI = stimCI(numCorrect);
  correctKernel = sum(cell2mat(U.hitKernel) .* double(U.stimCorrects), 1) / numCorrect;
  doOneKernelPlot(4, correctKernel, 'stim', plotStartMS, plotEndMS, plotTitle, ylabel, offset + correctCI, offset - correctCI);

  % fail kernel
  numFails = sum(U.stimFails);
  plotTitle = sprintf('Miss Kernel (n=%d)', numFails);
  missCI = stimCI(numFails);
  failKernel = sum(cell2mat(U.failKernel) .* double(U.stimFails), 1) / numFails;
  doOneKernelPlot(5, failKernel, 'stim', plotStartMS, plotEndMS, plotTitle, '', offset + missCI, offset - correctCI);

  % total kernel trial-weighted across all trials. We need to multiple the weighted sum by 2 because it is effectively
  % a mean of the hit and miss kernels, not a difference. By taking the mean, we lose the doubling that we should get
  % from the opposing effects.  This has been validated in simulations. 
  plotTitle = sprintf('Neuro-Behav Kernel (n=%d)', numCorrect + numFails);
  kernel = (sum((cell2mat(U.hitKernel) - 0.5) .* double(U.stimCorrects), 1) - ...
      sum((cell2mat(U.failKernel) - 0.5) .* double(U.stimFails), 1)) / sum(U.stimCorrects + U.stimFails) * 2;
  totalCI = stimCI(sum(U.stimCorrects + U.stimFails));
  kernel = doOneKernelPlot(6, kernel, 'stim', plotStartMS, plotEndMS, plotTitle, '', totalCI, -totalCI);
  sigmas = min(kernel) / -totalCI * 1.96;                   % use 95% CI to get SEMs of kernel peak
  
  % RT aligned kernel
  plotTitle = sprintf('RT Aligned (n=%d)', numCorrect);
  doOneKernelPlot(7, mean(cell2mat(U.RTKernel), 1), 'RT', plotRTStartMS, plotRTStartMS + plotEndMS - plotStartMS, ...
      plotTitle, ylabel, offset + correctCI, offset - correctCI);

  % Stim-RT aligned kernel
  plotTitle = sprintf('Stim-RT Aligned (n=%d)', numCorrect);
  doOneKernelPlot(10, mean(cell2mat(U.SRTKernel), 1), 'stimRT', plotRTStartMS, plotRTStartMS + plotEndMS - plotStartMS, ...
      plotTitle, ylabel, offset + correctCI, offset - correctCI);

  % Early kernel. There might be some sessions with no early releases, so we clear them out
  kernelLength = zeros(1, size(U, 1));
  for index = 1:size(U, 1)
      kernelLength(index) = size(U.earlyKernel{index}, 2);
  end
  validIndices = find(kernelLength > 1);
  if ~isempty(validIndices)
      earlyKernels = U.earlyKernel(validIndices);
      numEarlies = sum(U.stimEarlies(validIndices));
      plotTitle = sprintf('Early aligned (n=%d)', numEarlies);
      earlyCI = stimCI(numEarlies);
      doOneKernelPlot(8, mean(cell2mat(earlyKernels), 1), 'early', plotStartMS, plotEndMS, plotTitle, '', ...
              offset + earlyCI, offset - earlyCI);
  end

  % Random Kernel
  if numCorrect > 0 && numFails > 0
      plotTitle = sprintf('Random Kernel (n=%d)', numCorrect + numFails);
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
  earlyRTs = cat(2, U.earlyRTs{:});
  failRTs = cat(2, U.failRTs{:});
  doRTHistogramPlot(correctRTs, earlyRTs, failRTs, minRespTimeMS, maxRespTimeMS);
  doRTPDFPlot(correctRTs, earlyRTs, failRTs, minRespTimeMS, maxRespTimeMS)

  % Coordinate the scaling across plots
  sameYAxisScaling(4, 3, [4, 5, 7, 8, 10]);
  sameYAxisScaling(4, 3, [6, 9]);
   
  % display header text
  doHeader(U, limits);
end