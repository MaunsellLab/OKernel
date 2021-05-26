function makeMethodFigure()

  [U, dataDir, limits] = getSessionTable('Example');
  folderName = strcat(dataDir, ' Analysis/Figures/Figure XX Methods');
  if ~exist(folderName, 'dir')
   mkdir(folderName);
  end
  load(strcat(dataDir, U.animal, '/MatFiles/', U.date, '.mat'), 'file', 'trials');
  h = figure(10);
  set(h, 'Units', 'inches', 'Position', [25, 14.5, 8.0, 11.0]);
  clf;
  getStimProfiles(file, trials);
	saveas(gcf, strcat(folderName, '/OptoStim.pdf'));
  bootstraps = getCaseBootstraps(U, dataDir, limits.animal{1}, limits, true);
  limits.numBoot = 25;
  limits.aniNum = 1;
  limits.animal = '902';
  doOneBootFigure(U, dataDir, 'Methods', limits, bootstraps);
	saveas(gcf, strcat(folderName, '/Kernels.pdf'));

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
  meanPowers = [trials(:).meanPowerMW];                        % get power applied for each trial 
  row.meanPowerMW = mean(meanPowers);
  row.maxPowerMW = max(meanPowers);
  trialStructs = [trials(:).trial];
  eotCodes = zeros(1, length(trials));                        % one and only one RT and endTrial per trial
  for t = 1:length(trials)                               
      if length(trials(t).reactTimeMS) > 1
          trials(t).reactTimeMS = trials(t).reactTimeMS(1);   
      end
      if ~isfield(trials(t), 'trialEnd') || isempty(trials(t).trialEnd)
          trials(t).trialEnd = -1;
      end
      if length(trials(t).trialEnd) > 1
          trials(t).trialEnd = trials(t).trialEnd(1);   
      end
      eotCodes(t) = trials(t).trialEnd;
  end
  RTs = [trials(:).reactTimeMS];                              % get all trial RTs  
	preStimMS = [trialStructs(:).preStimMS];                  	% get preStim times for each trial  

  % calculate the overall d'
  theIndices = allIndices(trials, eotCodes, stimIndices);
  if sum(theIndices.correct | theIndices.fail | theIndices.early) < 10
    return;
  end
  % find the response interval and get a modified set of indices that limits hits to only that interval
	[respLimitsMS, theIndices, ~, ~] = getResponseLimits(file, trials, theIndices);
  row.RTWindowMS = diff(respLimitsMS);
  row.noStimCorrects = sum(theIndices.correct & ~stimIndices);
  row.noStimFails = sum(theIndices.fail & ~stimIndices);
  row.noStimEarlies = sum(theIndices.early & ~stimIndices); 
  row.numNoStim = row.noStimCorrects + row.noStimFails + row.noStimEarlies;  % don't count kEOTIgnored
  row.stimCorrects = sum(theIndices.correct & stimIndices);
  row.stimFails = sum(theIndices.fail & stimIndices);
  row.stimEarlies = sum(theIndices.early & stimIndices);
  row.numStim = row.stimCorrects + row.stimFails + row.stimEarlies;                  % don't count kEOTIgnored

  % find the performance across stim and nostim trials combined
  hitRate = sum(theIndices.correct) / (sum(theIndices.correct) + sum(theIndices.fail));
  rateEarly = earlyRate(file, trials, theIndices.correct, theIndices.fail, theIndices.early);
  row.pFA = 1.0 - exp(-rateEarly * row.RTWindowMS / 1000.0);
  row.pHit = (hitRate - row.pFA) / (1.0 - row.pFA);
  [row.dPrime, row.c] = dprime(row.pHit, row.pFA, true);
  
  % calculate the nostim trial d' using the all-trial pFA and the all-trial response window
  indices.correct = theIndices.correct & ~stimIndices;
  indices.fail = theIndices.fail & ~stimIndices;
  indices.early = theIndices.early & ~stimIndices;
  if sum(indices.correct | indices.fail | indices.early) > 0
    hitRate = sum(indices.correct) / (sum(indices.correct) + sum(indices.fail));
    rateEarly = earlyRate(file, trials, indices.correct, indices.fail, indices.early);
    row.noStimPFA = 1.0 - exp(-rateEarly * row.RTWindowMS / 1000.0);
    row.noStimPHit = (hitRate - row.pFA) / (1.0 - row.pFA);         % using overall pFA
    [row.noStimDPrime, row.noStimC] = dprime(row.noStimPHit, row.pFA, true);
  end
    
  % calculate the stim trial d' using the all-trial pFA and the all-trial response window
  indices.correct = theIndices.correct & stimIndices;
  indices.fail = theIndices.fail & stimIndices;
  indices.early = theIndices.early & stimIndices;
  if sum(indices.correct | indices.fail | indices.early) > 0
    hitRate = sum(indices.correct) / (sum(indices.correct) + sum(indices.fail));
    rateEarly = earlyRate(file, trials, indices.correct, indices.fail, indices.early);
    row.stimPFA = 1.0 - exp(-rateEarly * row.RTWindowMS / 1000.0);
    row.stimPHit = (hitRate - row.pFA) / (1.0 - row.pFA);   % using overall pFA
    [row.stimDPrime, row.stimC] = dprime(row.stimPHit, row.pFA, true);
  end
  
  % get the various kernels.  We use the indices set up in the previous block that include only stimulated trials.
  
%   [plotStartMS, plotEndMS, plotRTStartMS] = plotLimits();     % get the limits for the plots we will display  
  % get the hit kernel.  We use the indices set up to have only stimulated trials within the detected response window.
%   profiles = getTrialStimProfiles(trials(indices.correct), plotStartMS, plotEndMS, true, false);

  tIndices = find(indices.correct | indices.fail);
  
% The following will display all the individual trials, with a trial index.  Those index values can be used in the
% following alternative line to select files to plot in the figure.

%   profiles = screenProfiles(file, trials(indices.correct | indices.fail));
  selectedTrials = [11, 7, 5, 17, 8, 10];
  trialLabels = [sum(indices.correct | indices.fail), length(selectedTrials) - 1:-1:1];
  plotProfiles(file, trials(tIndices(selectedTrials)), trialLabels);
  selectedTrials = [11, 17, 10];
  trialLabels = [sum(indices.correct), length(selectedTrials) - 1:-1:1];
  plotProfileType(2, trials(tIndices(selectedTrials)), trialLabels);
  selectedTrials = [7, 5, 8];
  trialLabels = [sum(indices.fail), length(selectedTrials) - 1:-1:1];
  plotProfileType(3, trials(tIndices(selectedTrials)), trialLabels);
  
%   hitSums(:) = normSum(profiles);
%   hitSums(:) = sum(profiles, 1);
%   hitKernel = hitSums / row.stimCorrects;
%   row.hitKernel = {hitKernel};
%   hitCI = stimCI(row.stimCorrects);
% 
%   % get the RT aligned kernel
%   profiles = getTrialStimProfiles(trials(indices.correct), plotRTStartMS, plotRTStartMS + plotEndMS - plotStartMS, true, true);
% %   RTSum(:) = normSum(profiles);
%   RTSum(:) = sum(profiles, 1);
%   row.RTKernel = {RTSum / row.stimCorrects};
% 
%   % get the Stim & RT aligned hit kernel
%   profiles = getStimRTProfiles(trials(indices.correct), plotStartMS, plotEndMS);
% %   SRTSum(:) = normSum(profiles);
%   SRTSum(:) = sum(profiles, 1);
%   row.SRTKernel = {SRTSum / row.stimCorrects};
% 
%   % get the fail kernel
%   profiles = getTrialStimProfiles(trials(indices.fail), plotStartMS, plotEndMS, true, false);
%   failSum(:) = sum(profiles, 1);
%   failKernel = failSum / row.stimFails;
%   row.failKernel = {failKernel};
%   failCI = stimCI(row.stimFails);
% 
%   % Get the early kernel.  Eliminate trials that start before the end of the ramping stimulus.
% 	earlyIndices = indices.early & (preStimMS + RTs + plotRTStartMS > 200);
%   if row.stimEarlies > 0
%       profiles = getTrialStimProfiles(trials(earlyIndices), plotRTStartMS, plotRTStartMS + plotEndMS - plotStartMS, true, true);
%       earlySum(:) = normSum(profiles);
%       row.stimEarlies = size(profiles, 1);                  % getStimProfiles might reject some trials as too short
%       row.earlyKernel = {earlySum / row.stimEarlies};
%   end
% 
%   % If we have corrects and fails, save the hit, fail, total and random kernels.
%   row.kernelCI = sqrt(hitCI^2 + failCI^2);
%   if ~isempty(hitKernel) && ~isempty(failKernel)
%     row.peakMinMS = plotStartMS;
%     row.peakMaxMS = plotEndMS;
%     startIndex = row.peakMinMS - plotStartMS + 1;
%     endIndex = row.peakMaxMS - plotStartMS;
%     row.kernelPeak =  max(abs((hitKernel(startIndex:endIndex) - failKernel(startIndex:endIndex))) / row.kernelCI);
%     row.randomKernel = getRandomKernel(row.stimCorrects, row.stimFails, trialStructs(1).pulseDurMS, plotEndMS - plotStartMS);
%   end
% 
%   % add to the RT distributions
%   row.correctRTs = {[trials(theIndices.correct).reactTimeMS]};
%   row.earlyRTs = {[trials(theIndices.early).reactTimeMS]};
%   failRTs = [trials(theIndices.fail).reactTimeMS];
%   failRTs(failRTs < 0 | failRTs > 100000) = 100000;     % include fails in count, but don't let them display on plot
%   row.failRTs = {failRTs};
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
