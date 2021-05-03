function bootKernel
% Produce the kernel-STA convolutions with CIs.  Bootstrapped STA and kernel profiles are created, and then
% used to bootstrap convolutions CIs

% Get the bootstrap samples for the STA

  numBS = 10000;
  typeNames = {'Pyr', 'PV'};

  for t = 1:length(typeNames)
    typeName = typeNames{t};
    clear unitSTA;
    load(['STAs', typeName, 'Units'], 'unitSTA');
    fprintf(sprintf('Collecting %d bootstrap samples for %s STAs\n', numBS, typeName));
    profiles = bootstrp(numBS, @mean, unitSTA);
    bootSTA = zeros(numBS, 800);
    bootSTA(:, 1:500) = profiles - 0.5;
    PCs = prctile(bootSTA, [5, 50, 95]);
    CI05 = PCs(1,:);
    CI50 = PCs(2,:);
    CI95 = PCs(3,:);

    h = figure(1);
    set(h, 'Units', 'inches', 'Position', [25, 14.5, 8.0, 11.0]);
    clf;
    clf;
    subplot(3, 2, 1);
    plot(CI50);
    hold on;
    x = 1:size(CI50, 2);
    x2 = [x, fliplr(x)];
    fillCI = [CI05, fliplr(CI95)];
    fill(x2, fillCI, 'b', 'lineStyle', '-', 'edgeColor', 'b', 'edgeAlpha', 0.5, 'faceAlpha', 0.10);
    ax = axis;
    plot([ax(1), ax(2)], [0.0, 0.0], 'k--');
    plot([400, 400], [ax(3), ax(4)], 'k--');
    title(sprintf([typeName, ' STA weighted by units\n(%d units, %d bootstrap samples)'], size(unitSTA, 1), numBS));

    % Get the bootstrap samples for the kernel(s)

    dataDirName = '/Users/Shared/Data/OKernel/';
    tableDataName = [dataDirName ' Analysis/Processed Files.mat'];
    limits.minSessions = 8;   % at least 8 sessions for each animal
    limits.criterion = 0.0;   % no kernel significance criterion
    limits.animal = 'All';    % all animals
    limits.minTrials = 0;     % no trial minimum
    limits.minDec = 0;        % no negative change in behavioral performance
    limits.rampMS = 0;
    limits.stimStr = 'Steps';
    fprintf(sprintf('Processing %d bootstrap samples for %s rampMS %d\n', numBS, typeName, limits.rampMS));
    convOneStim(dataDirName, tableDataName, limits, bootSTA, typeName, 3);

    limits.rampMS = 500;
    limits.stimStr = 'Ramps';
    fprintf(sprintf('Processing %d bootstrap samples for %s rampMS %d\n', numBS, typeName, limits.rampMS));
    convOneStim(dataDirName, tableDataName, limits, bootSTA, typeName, 4);

    saveas(gcf, ['/Users/Shared/Data/OKernel/ Analysis/Figures/', 'BootConv ', typeName, limits.animal, '.pdf']);
    save(['STAs', typeName, 'Units'], 'unitSTA', 'CI50', 'CI05', 'CI95');

  end
end

%%
function convOneStim(dataDirName, tableDataName, limits, bootSTA, typeName, plotNum)

  numBS = size(bootSTA, 1);
  fileName = sprintf('profiles%03d.mat', limits.rampMS);  % load the stim profiles for the kernel
  if isfile(fileName)
    load(fileName, 'fileProfiles');
  else
    % get a table with all the valid sessions.
    [U, ~] = getSubset('normal', dataDirName, tableDataName, limits);
    fileProfiles = getOptoStim(U, dataDirName, tableDataName, limits);
    save(fileName, 'fileProfiles');
  end
  bootKernel = bootstrp(numBS, @mean, fileProfiles);      % bootstrap the kernel CIs     
  PCs = prctile(bootKernel, [5, 95]);                     % get the raw CIs
  
	% low-pass filtering was worked out using testFiltering.m
  sampleFreqHz = 1000;
  filterLP = designfilt('lowpassfir', 'PassbandFrequency', 90 / sampleFreqHz, ...
    'StopbandFrequency', 2 * 90 / sampleFreqHz, 'PassbandRipple', 1, 'StopbandAttenuation', 60, ...
    'DesignMethod','equiripple');

  subplot(3, 2, plotNum)
  profileMean = mean(fileProfiles);
  x = 1:length(profileMean);
  plot(x, filtfilt(filterLP, profileMean), 'b');
  hold on;
  CI05 = filtfilt(filterLP, PCs(1,:));
  CI95 = filtfilt(filterLP, PCs(2,:));
  x2 = [x, fliplr(x)];
  fillCI = [CI05, fliplr(CI95)];
  h = fill(x2, fillCI, 'b', 'lineStyle', '-', 'edgeColor', 'b', 'edgeAlpha', 0.5, 'faceAlpha', 0.10);
  set(h, 'faceAlpha', 0.10);
  ax = axis;
  plot([ax(1), ax(2)], [0, 0], 'k--');
  plot([1, 1] * mean(ax(1:2)), [ax(3), ax(4)], 'k--');
  title(sprintf('%s %s weighted by trial\n(%d trials, %d bootstrap samples)', ...
    typeName, limits.stimStr, size(fileProfiles, 1), numBS));
  
  bootConvs = zeros(numBS, 800);
  for i = 1:numBS
    c = xcorr(bootKernel(i,:), bootSTA(i,:));     % do row-wise convolutions for STA X kernel
    bootConvs(i,:) = c(401:1200);
  end
  PCs = prctile(bootConvs, [5, 50, 95]);          % get CIs on convolutions
  CI05 = PCs(1,:);
  CI50 = PCs(2,:);
  CI95 = PCs(3,:);
  subplot(3, 2, plotNum + 2);
  plot(CI50);
  hold on;
  x = 1:size(CI50, 2);
  x2 = [x, fliplr(x)];
  fillCI = [CI05, fliplr(CI95)];
  fill(x2, fillCI, 'b', 'lineStyle', '-', 'edgeColor', 'b', 'edgeAlpha', 0.5, 'faceAlpha', 0.10);
  ax = axis;
  plot([ax(1), ax(2)], [0.0, 0.0], 'k--');
  plot([400, 400], [ax(3), ax(4)], 'k--');
  title(sprintf('%s STA-%s Kernel Convolution', typeName, limits.stimStr));
end