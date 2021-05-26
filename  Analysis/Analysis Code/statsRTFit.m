function statsRTFit
%
% Get the r^2 for the logistic fits to the RT distributions for steps and ramps and combined
%
  [stepRespLimitsMS, stepRespLimitsR2] = stimRespStats('All Steps');
  [rampRespLimitsMS, rampRespLimitsR2] = stimRespStats('All Ramps');
  totalRespLimitsR2 = [stepRespLimitsR2', rampRespLimitsR2']';
  totalRespLimitsMS = [stepRespLimitsMS', rampRespLimitsMS']';
	showStats('Step and Ramp', totalRespLimitsR2, totalRespLimitsMS);
  [p, ~, stats] = ranksum(diff(stepRespLimitsMS'), diff(rampRespLimitsMS'));
  fprintf('Wilcoxon Rank Sum Test on Resp Window (Step v Ramp):\n  p = %.3f, z value %.2f\n', p, stats.zval);
end

function [respLimitsMS, respLimitsR2] = stimRespStats(subset)
  [U, dataDirName] = getSessionTable(subset);
  respLimitsMS = zeros(height(U), 2);
  respLimitsR2 = zeros(height(U), 1);
  for s = 1:height(U)
    if ~mod(s, 50)
      fprintf('Session %d of %d\n', s, height(U));
    end
    % Load and condition the data
    load(strcat(dataDirName, U.animal(s), '/MatFiles/', U.date(s)), 'file', 'trials');
    [stimIndices, trials] = getStimIndices(trials);
    if sum(stimIndices) == 0
      return
    end

    eotCodes = zeros(1, length(trials));                        % one and only one RT and endTrial per trial
    for t = 1:length(trials)                               
        if length(trials(t).reactTimeMS) > 1
            trials(t).reactTimeMS = trials(t).reactTimeMS(1);   
        end
        if ~isfield(trials(t), 'trialEnd') || isempty(trials(t).trialEnd)
            trials(t).trialEnd = -1;
        end
        eotCodes(t) = trials(t).trialEnd;
    end
    indices = allIndices(trials, eotCodes, stimIndices);
    [respLimitsMS(s, :), ~, ~, ~, fitStats] = getResponseLimits(file, trials, indices);
    respLimitsR2(s) = fitStats.rsquare;
  end
  showStats(subset, respLimitsR2, respLimitsMS);
end

function showStats(subset, limitsR2, limitsMS)
  fprintf('Response Window Statistics (%s) (n = %d)\n', subset, length(limitsR2));
  fprintf('  Average r^2 of logistic fit %.2f\n', mean(limitsR2));
%   fprintf('  Mean response window %.0f ms (SD %.0f ms)\n', mean(diff(limitsMS')), std(diff(limitsMS')));
  fprintf('  Median response window %.0f ms (IQR %.0f - %.0f ms)\n', prctile(diff(limitsMS'), [50, 25, 75]));
end