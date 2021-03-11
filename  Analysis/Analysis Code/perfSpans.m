function perfSpans
  % perfSpans --
  % examine how changing the number of contiguous trials affects the difference in hit rate between stim & no-stim
  % trials.
  
  minTrials = 10;
  dataDirName = '/Users/Shared/Data/OKernel/';
  tableDataName = [dataDirName ' Analysis/Processed Files.mat'];
  load(tableDataName, 'T');
  validRamp = (T.rampMS == 0 | T.rampMS == 500) & T.kernelCI > 0;   % empty entries have zero for kernelCI
  enoughHits = T.hits > minTrials;
  enoughMisses = T.misses > minTrials;
  validRow = validRamp & enoughHits & enoughMisses;
  U = T(validRow, :);
  U = withoutTrainingDays(U);
  doPlots(U, dataDirName, minTrials);
end

%%
function doPlots(U, dataDirName, minTrials)

  display(U);
  h = figure(1);
  set(h, 'Units', 'inches', 'Position', [25, 1.25, 8.5, 11]);
  for f = 1:height(U)
    fprintf('%s %s\n', U.animal(f), U.date(f));
    if mod(f - 1, 12) == 0
      if f > 1
        saveas(h, sprintf('~/Desktop/Figures/Figure %02.0f.pdf', f / 12));
      end
      clf;
    end
    subplot(4, 3, mod(f - 1,12) + 1);
    clear trials;
    load(sprintf('%s%s%s%s', dataDirName, U.animal(f), '/MatFiles/', U.date(f)), 'trials');
    for t = 1:length(trials)                                    % one, only one RT and endTrail per trial
      if ~isfield(trials(t), 'trialEnd') || isempty(trials(t).trialEnd)
        trials(t).trialEnd = -1;
      end
      if ~isfield(trials(t), 'meanPowerMW') || isempty(trials(t).meanPowerMW)
        trials(t).meanPowerMW = -1;
      end
    end
    % get the overall unstimulated hit rate during the stimulated portion of the session
    meanPower = [trials(:).meanPowerMW];                        % get power applied on each trial
    firstStimIndex = find(meanPower > 0, 1);                    % first stimulated trial
    lastStimIndex = find(meanPower > 0, 1, 'last');             % last stimulated trial
    trials = trials(firstStimIndex:lastStimIndex);              % take the stimulated portion only
    meanPower = [trials(:).meanPowerMW];                        % get power applied on each trial
    eotCodes = [trials(:).trialEnd];                            % get eotCodes
    hitTrials = sum(meanPower == 0 & eotCodes == 0);
    missTrials = sum(meanPower == 0 & eotCodes == 2);
    unstimHitRate = hitTrials / (hitTrials + missTrials);
    trials = trials(meanPower > 0 & (eotCodes == 0 | eotCodes == 2));   % take the stimulated trials only
    numTrials = length(trials);
    if numTrials < minTrials
      continue;
    end
    eotCodes = [trials(:).trialEnd];                            % get eotCodes
    hitTrials = eotCodes == 0;
    missTrials = eotCodes == 2;
    spanRates = ones(1, numTrials);                             % decrement in hit rate
    simRates = ones(1, numTrials);                              % simulated in hit rate
    pHit = sum(hitTrials) / (sum(hitTrials) + sum(missTrials));
    for spanLength = numTrials:-1:1                             % for all spans between max and minTrials
      for s = 1:numTrials - spanLength + 1                      % for all possible spans of this length
        hits = sum(hitTrials(s:s + spanLength - 1));
        misses = sum(missTrials(s:s + spanLength - 1));
        hitRate = hits / (hits + misses);
%         fprintf('spanLength %2d, s %2d, hits %2d, misses %2d, rate %.2f\n', spanLength, s, hits, misses, hitRate);
        spanRates(spanLength) = min(spanRates(spanLength), hitRate);
        simRates(spanLength) = min(simRates(spanLength), sum(rand(1, spanLength) < pHit) / spanLength);
      end
    end
    plot(1:numTrials, spanRates(1:end), 'b');
    hold on;
    plot(1:numTrials, simRates(1:end), 'b:');
    plot([minTrials, length(trials)], [unstimHitRate, unstimHitRate], 'r');
    plot([minTrials, length(trials)], [unstimHitRate - 0.15, unstimHitRate - 0.15], 'r:');
    ylim([0, 1]);
    title(sprintf('%s %s', U.animal(f), U.date(f)));
    if (f >= 10)
      xlabel('Number of Trials Included');
    end
    if mod(f-1, 3) == 0
      ylabel('Hit Rate');
    end
    drawnow;
  end
end

%%
function U = withoutTrainingDays(U)

	trainingDays = ...
    (U.animal == '902' & U.date < '2019-09-13') | ...
    (U.animal == '905' & U.date < '2019-09-24') | ...
    (U.animal == '1112' & U.date < '2020-02-04') | ...
    (U.animal == '1145' & U.date < '2020-02-10') | ...
    (U.animal == '1150' & U.date < '2020-01-16') | ...
    (U.animal == '1218' & U.date < '2020-03-23') | ...
    (U.animal == '1220' & U.date < '2020-04-07') | ...
    (U.animal == '1223' & U.date < '2020-02-16') | ...
    (U.animal == '1257' & U.date < '2020-04-04');    
  U = U(~trainingDays, :);

end
