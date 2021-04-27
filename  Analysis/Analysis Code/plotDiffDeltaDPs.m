function plotDiffDeltaDPs
% Plot different delta-d' subsets for the data

%{
Try taking each animal and splitting them into high/low delta-d' sessions (median split).
Might also do this by trials instead, to keep things balanced on significance. 
%}
	dataDirName = '/Users/Shared/Data/OKernel/';
  load([dataDirName ' Analysis/Mat Files/masterTable.mat'], 'T');

  rampMS = 0;
  animals = {'902', '1112', '1145', '905', '1223'};
  limits = setLimits('All Steps');
  limits.rampMS = rampMS;
  limits.animals = animals;
  limits.numBoot = 10;
  T = selectUsingLimits(T, limits);

% For these analyses, we need valid d' measures.
  T = T(~isnan(T.dPrime) & ~isnan(T.stimDPrime) & ~isnan(T.noStimDPrime), :);

%   rampMS = 500;
%   animals = {'902', '1112', '1150'};
  
  h = figure(2);
  set(h, 'Units', 'inches', 'Position', [25, 10, 8.5, 11]);
  clf;  
  plotKernelVDeltaDP(T);
 
  h = figure(1);
  set(h, 'Units', 'inches', 'Position', [25, 1.25, 8.5, 11]);
  clf;  
  [smallDPs, bigDPs] = splitSessions(T, animals);
	doFigure(dataDirName, smallDPs, bigDPs, limits);
	saveas(gcf, [dataDirName, ' Analysis/Figures/Big-Small DPs.pdf']);
end

%%
function doFigure(dataDirName, smallDPs, bigDPs, limits)
  % Compile and plot the kernels

 	smallBootstraps = getCaseBootstraps(smallDPs, dataDirName, 'Small Delta d''', limits);
	bigBootstraps = getCaseBootstraps(bigDPs, dataDirName, 'Large Delta d''', limits);
  [plotStartMS, plotEndMS, ~] = plotLimits();
  ylabel = 'Normalized Power';
  limits.yAxis = 0.5;
  
  % hit kernel
  numHits = size(smallBootstraps.hitProfiles, 1);
  plotTitle = sprintf('Hit Kernel (n=%d)', numHits);
  subplot(4, 3, 4);
  smallHitCIs = doOneBootPlot(smallBootstraps.hitProfiles / 2 + 0.5, limits, 'stim', plotStartMS, plotEndMS, plotTitle, ylabel);
%   save([dataDirName, ' Analysis/Mat Files/', dataName, ' ', limits.animal{1}, ' Hit Kernel'], 'CIs');
  
  % miss kernel
  numMisses = size(smallBootstraps.missProfiles, 1);
  plotTitle = sprintf('Miss Kernel (n=%d)', numMisses);
  subplot(4, 3, 5);
  smallMissCIs = doOneBootPlot(smallBootstraps.missProfiles / 2 + 0.5, limits, 'stim', plotStartMS, plotEndMS, plotTitle, '');

  % full kernel
  plotTitle = sprintf('Full Kernel (n=%d)', numHits + numMisses);
  limits.yAxis = 0.0;
  hitMissBoot = [smallBootstraps.hitProfiles; -smallBootstraps.missProfiles];
  subplot(4, 3, 6);
  smallFullCIs = doOneBootPlot(hitMissBoot, limits, 'stim', plotStartMS, plotEndMS, plotTitle, '');
  headerPlot(1, smallDPs, limits);
  
  % hit kernel
  numHits = size(bigBootstraps.hitProfiles, 1);
  plotTitle = sprintf('Hit Kernel (n=%d)', numHits);
  subplot(4, 3, 7);
  bigHitCIs = doOneBootPlot(bigBootstraps.hitProfiles / 2 + 0.5, limits, 'stim', plotStartMS, plotEndMS, plotTitle, ylabel);
%   save([dataDirName, ' Analysis/Mat Files/', dataName, ' ', limits.animal{1}, ' Hit Kernel'], 'CIs');
  
  % miss kernel
  numMisses = size(bigBootstraps.missProfiles, 1);
  plotTitle = sprintf('Miss Kernel (n=%d)', numMisses);
  subplot(4, 3, 8);
  bigMissCIs = doOneBootPlot(bigBootstraps.missProfiles / 2 + 0.5, limits, 'stim', plotStartMS, plotEndMS, plotTitle, '');

  % full kernel
  plotTitle = sprintf('Full Kernel (n=%d)', numHits + numMisses);
  limits.yAxis = 0.0;
  hitMissBoot = [bigBootstraps.hitProfiles; -bigBootstraps.missProfiles];
  subplot(4, 3, 9);
  bigFullCIs = doOneBootPlot(hitMissBoot, limits, 'stim', plotStartMS, plotEndMS, plotTitle, '');
  headerPlot(2, bigDPs, limits);

  pairedBootPlots(10, smallHitCIs, bigHitCIs, limits, 'Hit Kernels', '');
  pairedBootPlots(11, smallMissCIs, bigMissCIs, limits, 'Miss Kernels', '');
  pairedBootPlots(12, smallFullCIs, bigFullCIs, limits, 'Full Kernels', '');

  % Coordinate the scaling across plots
  sameYAxisScaling(4, 3, [4, 5, 7, 8, 10, 11], [0.40, 0.60]);
  sameYAxisScaling(4, 3, [6, 9, 12], [-0.15, 0.05]);

end

%%
function headerPlot(plotIndex, U, limits)
%
% Display the header information for an OKernel analysis page
%
  headerText = cell(1, 1);
  if limits.rampMS == 0
    headerText{1} = sprintf('Visual Contrast Step');
  else
    headerText{1} = sprintf('Visual Contrast Ramp %d ms', limits.rampMS);
  end
  if strcmp(limits.animal, 'All')
    headerText{end + 1} = sprintf('%d sessions from %d animals', size(U, 1), length(unique(U.animal)));
    animalList = [];
    animals = unique(U.animal);
    for a = 1:length(animals)
      animalList = strcat(animalList, sprintf('  %s (%d),', animals{a}, sum(U.animal == animals{a})));
      if ~mod(a, 3)
        headerText{end + 1} = animalList; %#ok<AGROW>
        animalList = [];
      end
    end
    if mod(a, 4)
      headerText{end + 1} = animalList; 
    end
  else
    headerText{end + 1} = sprintf('%d sessions from Animal %s', size(U, 1), limits.animal);
  end
  if limits.minDeltaDPrime == -1
    headerText{end + 1} = 'No required delta-d'' with opto';
  else
    headerText{end + 1} = sprintf('Delta d'' >=%.2f', limits.minDeltaDPrime);
  end
  headerText{end + 1} = sprintf('Avg (SD) d'' noStim: %.2f (%.2f)', nanmean(U.noStimDPrime), nanstd(U.noStimDPrime));
  headerText{end + 1} = sprintf('Avg (SD) d''   stim: %.2f (%.2f)', nanmean(U.stimDPrime), nanstd(U.stimDPrime));
  headerText{end + 1} = sprintf('Avg (SD) delta d'': %.2f (%.2f)', nanmean(U.noStimDPrime - U.stimDPrime), ...
    nanstd(U.noStimDPrime - U.stimDPrime));
  headerText{end + 1} = sprintf('Avg (SD) power (mW): %.3f (%.3f)', mean(U.meanPowerMW), std(U.meanPowerMW));
  axisHandle = subplot(4, 3, plotIndex);						% default axes are 0 to 1
  set(axisHandle, 'visible', 'off');
  if plotIndex == 1
    text(0.00, 1.25, 'Small Delta d''', 'FontWeight', 'bold', 'FontSize', 16);
  else
    text(0.00, 1.25, 'Large Delta d''', 'FontWeight', 'bold', 'FontSize', 16);
  end
  text(0.00, 1.10, headerText, 'VerticalAlignment', 'top');
end

%%
function pairedBootPlots(plotIndex, smallHitCIs, bigHitCIs, limits, plotTitle, yLabel)

  subplot(4, 3, plotIndex);
  cla reset;  
	x = 1:size(smallHitCIs, 2);
  plot(x, smallHitCIs(2, :), 'b');
  hold on;
  x2 = [x, fliplr(x)];
  fillCI = [smallHitCIs(1, :), fliplr(smallHitCIs(3, :))];
  fill(x2, fillCI, 'b', 'lineStyle', '-', 'edgeColor', 'b', 'edgeAlpha', 0.5, 'faceAlpha', 0.10);
  plot(x, bigHitCIs(2, :), 'r');
  hold on;
  x2 = [x, fliplr(x)];
  fillCI = [bigHitCIs(1, :), fliplr(bigHitCIs(3, :))];
  fill(x2, fillCI, 'r', 'lineStyle', '-', 'edgeColor', 'r', 'edgeAlpha', 0.5, 'faceAlpha', 0.10);

  bins = length(smallHitCIs);
  ax = gca;
  xlim(ax, [0, bins]);
  ax.XGrid = 'on';
  [plotStartMS, plotEndMS, ~] = plotLimits();
	plot([0, bins], [limits.yAxis, limits.yAxis], 'k-');
  set(gca,'XTick', [0, -plotStartMS, -plotStartMS + 100, bins]);
  set(gca, 'XTickLabel', {sprintf('%d', plotStartMS), '0', '', sprintf('%d', plotEndMS)});
  xlabel('Time Relative to Stimulus');
  ylabel(yLabel);
  title(plotTitle);
  hold off;
end

%%
function plotKernelVDeltaDP(T)

  subplot(2, 1, 2);
  kernelPeak = -T.kernelPeak .* T.kernelCI;
  deltaDPrime = T.stimDPrime - T.noStimDPrime;
  plot(deltaDPrime, kernelPeak, 'o');
  xlabel('Delta d''');
  ylabel('Kernel Peak (normalized power)'); 
  ylim([-2, 0]);
  
  subplot(2, 1, 1);
  meanPowerMW = T.meanPowerMW;
  plot(meanPowerMW, deltaDPrime, 'o');
  ylabel('Delta d''');
  xlabel('Mean Power (mW)'); 
  ylim([-2, 0]);
end

%%
function [smallDPs, bigDPs] = splitSessions(T, animals)

  smallDPs = []; bigDPs = [];
  for a = 1:length(animals)
    A = T(T.animal == animals{a}, :);         % extract rows for this animals
    A.dPrime = A.stimDPrime - A.noStimDPrime; % replace dPrime with delta dPrime
    A = sortrows(A, 'dPrime');
    targetSplit = sum(A.numStim) / 2;
    row = 1; bigStimNum = 0;
    while bigStimNum < targetSplit
      bigStimNum = bigStimNum + A.numStim(row);
      row = row + 1;
    end
    % Check whether we'd be closer to the split if we stopped one session earlier
    if bigStimNum - targetSplit > targetSplit - (bigStimNum - A.numStim(row - 1))
      row = row - 1;
    end
    bigDPs = [bigDPs; A(1:row - 1, :)]; %#ok<AGROW>
    smallDPs = [smallDPs; A(row:end, :)]; %#ok<AGROW>
  end
end
