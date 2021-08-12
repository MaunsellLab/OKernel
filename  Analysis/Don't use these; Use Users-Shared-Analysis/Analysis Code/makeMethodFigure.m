function makeMethodFigure()

  [U, dataDir, limits] = getSessionTable('Example');
  folderName = strcat(dataDir, ' Analysis/Figures/Figure XX Methods');
  if ~exist(folderName, 'dir')
   mkdir(folderName);
  end
  load(strcat(dataDir, U.animal, '/MatFiles/', U.date, '.mat'), 'file', 'trials');
  h = figure(10);
  set(h, 'Units', 'inches', 'Position', [25, 14.5, 8.0, 10.5]);
  clf;
  getStimProfiles(file, trials);
% 	saveas(gcf, strcat(folderName, '/OptoStim.pdf'));
  exportgraphics(gcf, strcat(folderName, '/OptoStim.pdf'),'contentType', 'vector');
  % Get the plots for this figure, which includes the hit, miss and full kernels
  bootstraps = getCaseBootstraps(U, dataDir, limits.animal{1}, limits, true);
  limits.numBoot = 100;
  limits.aniNum = 1;
  limits.animal = '902';
  doOneBootFigure(U, dataDir, 'Methods', limits, bootstraps);
  sameYAxisScaling(4, 3, [4, 5], [0, 1]);         % force the kernels to have the same vertical scaling
	sameYAxisScaling(4, 3, 6, [-0.5, 0.5]);
% 	saveas(gcf, strcat(folderName, '/Kernels.pdf'));
  exportgraphics(gcf, strcat(folderName, '/Kernels.pdf'),'contentType', 'vector');
end

%%
function getStimProfiles(file, trials)
%
% Compute kernels for one session
%
  [stimIndices, trials] = getStimIndices(trials);
  if sum(stimIndices) == 0
    return
  end  
  % get the raw indices for EOTs in the data file
  [theIndices, trials] = allIndices(trials);
  if sum(theIndices.correct | theIndices.fail | theIndices.early) < 10
    return;
  end
  % find the response interval and get a modified set of indices that limits hits to only that interval
	[~, theIndices, ~, ~] = getResponseLimits(file, trials, theIndices);
      
  % calculate the stim trial d' using the all-trial pFA and the all-trial response window
  indices.correct = theIndices.correct & stimIndices;
  indices.fail = theIndices.fail & stimIndices;
  
% The following will display all the individual trials, with a trial index.  Those index values can be used in the
% following alternative line to select files to plot in the figure.

%   profiles = screenProfiles(file, trials(indices.correct | indices.fail));
  
%% The following will display a selected subset of the trials, selected using the call to screenProfiles

  tIndices = find(indices.correct | indices.fail);
  selectedTrials = [11, 7, 5, 17, 8, 10];
  trialLabels = [sum(indices.correct | indices.fail), length(selectedTrials) - 1:-1:1];
  plotProfiles(file, trials(tIndices(selectedTrials)), trialLabels);
  selectedTrials = [11, 17, 10];
  trialLabels = [sum(indices.correct), length(selectedTrials) - 1:-1:1];
  plotProfileType(2, trials(tIndices(selectedTrials)), trialLabels);
  selectedTrials = [7, 5, 8];
  trialLabels = [sum(indices.fail), length(selectedTrials) - 1:-1:1];
  plotProfileType(3, trials(tIndices(selectedTrials)), trialLabels);
  
end

%%
function stimProfiles = screenProfiles(file, trials) %#ok<DEFNU>
% getTrialStimProfiles()
% Return a millisecond-resolution copy of the optogenetic stimulus profile from each of a set of trials. All trials
% are aligned with the preStimTime
  
  numTrials = size(trials, 2);
  if numTrials == 0
      stimProfiles = [];
      return
  end
  trialStructs = [trials(:).trial];
  preStimMS = [trialStructs(:).preStimMS];
  visualDurMS = [trialStructs(:).visualDurMS];
  endTimeMS = max(preStimMS) + max(visualDurMS);
  startTimeMS = 0;  
  stimProfiles = zeros(numTrials, endTimeMS);
  figure(10);
  for t = 1:numTrials
    meanPower = trials(t).meanPowerMW;
    if meanPower == 0
      stimProfiles(t, :) = 0.5 * ones(1, endTimeMS - startTimeMS);
      continue;
    end
    x = trials(t).optoStepTimesMS(10:end);      % don't plot the ramp
    reactTimeMS = trials(t).reactTimeMS - x(1);
    stimOnTimeMS = -x(1);
    rewardedLimitMS = file.rewardedLimitMS - x(1);
    x = x - x(1);
    p = trials(t).optoStepPowerMW(10:end);
    limitIndex = find(x > rewardedLimitMS, 1);
    x = x(1:limitIndex);
    p = p(1:limitIndex);
    offset = mod(t-1, 10) * trials(t).meanPowerMW * 2.5;
    if offset == 0
      clf;
    end
    hold on;
    for i = 1:length(x) - 1
      fill([x(i), x(i), x(i + 1), x(i + 1)], [0, p(i), p(i), 0] + offset, [0.8, 0.8, 1], 'lineStyle', 'none');
    end
    plot([x(1), x(1), x(2)], [0, p(1), p(1)] + offset, 'k');
    for i = 2:length(x) - 1
      plot([x(i), x(i), x(i + 1)], [p(i - 1), p(i), p(i)] + offset, 'k');
    end
    plot([x(i), x(i)], [p(i), 0] + offset, 'k');
    plot([stimOnTimeMS, stimOnTimeMS], [0, 2.1 * trials(t).meanPowerMW] + offset, 'b');
    plot([rewardedLimitMS, rewardedLimitMS], [0, 2.1 * trials(t).meanPowerMW] + offset, 'k');
    if trials(t).trialEnd == 0
      plot([reactTimeMS, reactTimeMS], [0, 2.1 * trials(t).meanPowerMW] + offset, 'r');
    end
    text(x(end) + 25, offset + trials(t).meanPowerMW, sprintf('%d - %d', t, trials(t).trialEnd));
  end
end

%%
function plotProfiles(file, trials, trialLabels)
% 
% Plot a series of trials showing the full trial duration for hits and miss trials
  
  subplot(3, 1, 1);
  numTrials = size(trials, 2);
  for t = 1:numTrials
    x = trials(t).optoStepTimesMS(10:end);                        % don't plot the ramp
    reactTimeMS = trials(t).reactTimeMS - x(1);
    stimOnTimeMS = -x(1);
    rewardedLimitMS = file.rewardedLimitMS - x(1);
    x = x - x(1);
    p = trials(t).optoStepPowerMW(10:end);
    limitIndex = find(x > rewardedLimitMS, 1);
    x = x(1:limitIndex);
    p = p(1:limitIndex);
    offsetFactor = 3;
    if mod(t-1, 10) == 0
      offset = 0.0 * offsetFactor + trials(t).meanPowerMW;
    else
      offset = (mod(t-1, 10) + 1) * trials(t).meanPowerMW * offsetFactor + trials(t).meanPowerMW;
    end
    hold on;
    for i = 1:length(x) - 1
      fill([x(i), x(i), x(i + 1), x(i + 1)], [0, p(i), p(i), 0] + offset, [0.8, 0.8, 0.8], 'lineStyle', 'none');
    end
    plot([x(1), x(1), x(2)], [0, p(1), p(1)] + offset, 'k');
    for i = 2:length(x) - 1
      plot([x(i), x(i), x(i + 1)], [p(i - 1), p(i), p(i)] + offset, 'k');
    end
    plot([x(i), x(i)], [p(i), 0] + offset, 'k');
    lineFactor = 2.5;
    plot([stimOnTimeMS, stimOnTimeMS], [0, lineFactor * trials(t).meanPowerMW] + offset, 'b', 'lineWidth', 2);
    plot([rewardedLimitMS, rewardedLimitMS], [0, lineFactor * trials(t).meanPowerMW] + offset, 'k', 'lineWidth', 2);
    ylim([0, (numTrials + 1.5) * trials(t).meanPowerMW * offsetFactor + trials(t).meanPowerMW]);
    xlim([-250, 2750]);
    if trials(t).trialEnd == 0
      plot([reactTimeMS, reactTimeMS], [0, lineFactor * trials(t).meanPowerMW] + offset, 'r', 'lineWidth', 2);
    end
    if trials(t).trialEnd == 0
      text(x(end) + 25, offset + trials(t).meanPowerMW, 'Hit');
    else
      text(x(end) + 25, offset + trials(t).meanPowerMW, 'Miss');
    end
    text(-100, offset + trials(t).meanPowerMW, sprintf('%d', trialLabels(t)), 'horizontalAlignment', 'right');
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
  end
end

%%
function plotProfileType(plotIndex, trials, trialLabels)
% 
% Make a plot of hit or miss trials (separately) at the same scale as the combined subplot.  Shows only +-400 ms
% around stimulus onset to match the behavioral kernels
  
  subplot(3, 1, plotIndex);
  cla;
  numTrials = size(trials, 2);
  for t = 1:numTrials
    firstIndex = find(trials(t).optoStepTimesMS > -400, 1) - 1;
    trials(t).optoStepTimesMS(firstIndex) = -400;
    lastIndex = find(trials(t).optoStepTimesMS > 400, 1) - 1;
    trials(t).optoStepTimesMS(lastIndex) = 400;
    trials(t).optoStepPowerMW(lastIndex) = trials(t).optoStepPowerMW(lastIndex - 1);
    x = trials(t).optoStepTimesMS(firstIndex:lastIndex);                        % don't plot the ramp
    reactTimeMS = trials(t).reactTimeMS;
    stimOnTimeMS = 0;
    p = trials(t).optoStepPowerMW(firstIndex:lastIndex);
    offsetFactor = 3;
    if mod(t-1, 10) == 0
      offset = 0.0 * offsetFactor + trials(t).meanPowerMW;
    else
      offset = (mod(t-1, 10) + 1) * trials(t).meanPowerMW * offsetFactor + trials(t).meanPowerMW;
    end
    hold on;
    for i = 1:length(x) - 1
      fill([x(i), x(i), x(i + 1), x(i + 1)], [0, p(i), p(i), 0] + offset, [0.8, 0.8, 0.8], 'lineStyle', 'none');
    end
    plot([x(1), x(1), x(2)], [0, p(1), p(1)] + offset, 'k');
    for i = 2:length(x) - 1
      plot([x(i), x(i), x(i + 1)], [p(i - 1), p(i), p(i)] + offset, 'k');
    end
    plot([x(i), x(i)], [p(i), 0] + offset, 'k');
    lineFactor = 2.5;
    plot([stimOnTimeMS, stimOnTimeMS], [0, lineFactor * trials(t).meanPowerMW] + offset, 'b', 'lineWidth', 2);
    if trials(t).trialEnd == 0
      plot([reactTimeMS, reactTimeMS], [0, lineFactor * trials(t).meanPowerMW] + offset, 'r', 'lineWidth', 2);
    end
    if trials(t).trialEnd == 0
      text(x(end) + 25, offset + trials(t).meanPowerMW, 'Hit');
    else
      text(x(end) + 25, offset + trials(t).meanPowerMW, 'Miss');
    end
    text(-500, offset + trials(t).meanPowerMW, sprintf('%d', trialLabels(t)), 'horizontalAlignment', 'right');
    set(gca,'xtick',[-400, 0, 400]);
    set(gca,'xticklabel',{'-400', '0', '400'});
    set(gca,'ytick',[]);
    set(gca,'yticklabel',[]);
    xlim([-650, 2350]);
    ylim([0, (numTrials + 4.5) * trials(t).meanPowerMW * offsetFactor + trials(t).meanPowerMW]);
  end
end
