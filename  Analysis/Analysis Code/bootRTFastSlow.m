function bootRTFastSlow

  numBS = 10000;
  h = figure(1);
  set(h, 'Units', 'inches', 'Position', [25, 14.5, 8.0, 11.0]);
  clf;

% Get the bootstrap samples for the kernel(s)
 
  dataDirName = '/Users/Shared/Data/OKernel/';
	tableDataName = [dataDirName ' Analysis/Processed Files.mat'];
  limits.minSessions = 8;   % at least 8 sessions for each animal
  limits.criterion = 0.0;   % no kernel significance criterion
  limits.animal = 'All';    % all animals
  limits.minTrials = 0;     % no trial minimum
  limits.minDec = 0;        % no negative change in behavioral performance
  rampMS = [0, 500];
  limits.stimStr = {'Steps', 'Ramps'};
  for r = 1:2
    limits.rampMS = rampMS(r);
    doOneStim(dataDirName, tableDataName, limits, numBS, r);
  end
  plots = 1:10;
  yLimits = sameYAxisScaling(5, 2, plots);
  for p = plots
    subplot(5, 2, p);
    hold on;
    plot([0, 800], [0, 0], 'k--');
    plot([400, 400], yLimits, 'k--');
  end
  saveas(gcf, ['/Users/Shared/Data/OKernel/ Analysis/Figures/', 'BootRT ', limits.animal, '.pdf']);
end

%%
function doOneStim(dataDirName, tableDataName, limits, numBS, rep)

	fprintf(sprintf('Processing %d bootstrap samples for rampMS %d\n', numBS, limits.rampMS));
  fileName = sprintf('profiles%03d.mat', limits.rampMS);
  if isfile(fileName)
    load(fileName, 'fileProfiles');
  else
    % get a table with all the valid sessions.
    [U, ~] = getSubset('normal', dataDirName, tableDataName, limits);
    fileProfiles = getOptoStim(U, dataDirName);
    save(fileName, 'fileProfiles');
  end
  bootKernel = bootstrp(numBS, @mean, fileProfiles);
  titleStr = sprintf('All trials with %s\n(%d trials, %d bootstrap samples)', ...
                            limits.stimStr{rep}, size(fileProfiles, 1), numBS);
  plotKernel(bootKernel, rep, titleStr);

  % Get the fast and slow hits (without misses)
  
	fprintf(' Processing fast and slow RTs\n');
  fileName = sprintf('profilesFastSlow%03d.mat', limits.rampMS);
  if isfile(fileName)
    load(fileName, 'fileProfilesFast', 'fileProfilesSlow');
  else
    % get a table with all the valid sessions.
    [U, ~] = getSubset('normal', dataDirName, tableDataName, limits);
    [fileProfilesFast, fileProfilesSlow] = getOptoStimFastSlow(U, dataDirName);
    save(fileName, 'fileProfilesFast', 'fileProfilesSlow');
  end
  
  bootKernelFast = bootstrp(numBS, @mean, fileProfilesFast);
  titleStr = sprintf('Fast RT trials with %s\n(%d trials, %d bootstrap samples)', ...
                            limits.stimStr{rep}, size(fileProfilesFast, 1), size(bootKernelFast, 1));
  plotKernel(bootKernelFast, rep + 6, titleStr);
  
  bootKernelSlow = bootstrp(numBS, @mean, fileProfilesSlow);
  titleStr = sprintf('Slow RT trials with %s\n(%d trials, %d bootstrap samples)', ...
                            limits.stimStr{rep}, size(fileProfilesSlow, 1), size(bootKernelSlow, 1));
  plotKernel(bootKernelSlow, rep + 8, titleStr);

  % Get the fast and slow hits (with misses)
  
  fprintf(' Processing fast and slow RTs with Misses\n');
  fileName = sprintf('profilesFastSlowMiss%03d.mat', limits.rampMS);
  if isfile(fileName)
    load(fileName, 'fileProfilesFastMiss', 'fileProfilesSlowMiss');
  else
    % get a table with all the valid sessions.
    [U, ~] = getSubset('normal', dataDirName, tableDataName, limits);
    [fileProfilesFastMiss, fileProfilesSlowMiss] = getOptoStimFastSlowMiss(U, dataDirName);
    save(fileName, 'fileProfilesFastMiss', 'fileProfilesSlowMiss');
  end
  
  bootKernelFastMiss = bootstrp(numBS, @mean, fileProfilesFastMiss);
  titleStr = sprintf('Fast RT & Miss trials with %s\n(%d trials, %d bootstrap samples)', ...
                            limits.stimStr{rep}, size(fileProfilesFastMiss, 1), size(bootKernelFastMiss, 1));
  plotKernel(bootKernelFastMiss, rep + 2, titleStr);
  
  bootKernelSlowMiss = bootstrp(numBS, @mean, fileProfilesSlowMiss);
  titleStr = sprintf('Slow RT trials with %s\n(%d trials, %d bootstrap samples)', ...
                            limits.stimStr{rep}, size(fileProfilesSlowMiss, 1), size(bootKernelSlowMiss, 1));
  plotKernel(bootKernelSlowMiss, rep + 4, titleStr);

end

function plotKernel(bootKernel, plotNum, titleStr)

	% low-pass filtering was worked out using testFiltering.m
  sampleFreqHz = 1000;
  filterLP = designfilt('lowpassfir', 'PassbandFrequency', 90 / sampleFreqHz, ...
    'StopbandFrequency', 2 * 90 / sampleFreqHz, 'PassbandRipple', 1, 'StopbandAttenuation', 60, ...
    'DesignMethod','equiripple');

  PCs = prctile(bootKernel, [5, 50, 95]);
  CI05 = filtfilt(filterLP, PCs(1,:));
  CI50 = filtfilt(filterLP, PCs(2,:));
  CI95 = filtfilt(filterLP, PCs(3,:));
  subplot(5, 2, plotNum)
  x = 1:length(CI50);
  plot(x, CI50, 'b');
  hold on;
  x2 = [x, fliplr(x)];
  fillCI = [CI05, fliplr(CI95)];
  h = fill(x2, fillCI, 'b', 'lineStyle', '-', 'edgeColor', 'b', 'edgeAlpha', 0.5, 'faceAlpha', 0.10);
  title(titleStr);
  xticks(0:100:800);
  xticklabels({'-400', '', '-200', '', '0', '', '200', '', '400'});
  set(gca, 'xgrid', 'on');
end