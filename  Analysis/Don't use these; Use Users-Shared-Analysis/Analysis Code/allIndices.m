function [indices, trials] = allIndices(trials, eotCodes)
%
% Return logical arrays specifying valid indices for hit, fail and early trials.  The lists are modified
% according to trialType, which can select all, stimulated, or unstimulated trials.

  if nargin < 2
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
  end
	indices.correct = eotCodes == 0;                              % get hit trials
  indices.early = eotCodes == 1;                                % get early trials
  indices.fail = eotCodes == 2;                                 % get fail trials
  
  % for earlies, we must eliminate those that don't have enough time before the FA to avoid sampling before 
  % the end of the ramp (which is fixed at 200 ms). 
  trialStructs = [trials(:).trial];                           % extract trialStructs to an array for access
  [~, ~, plotRTStartMS] = plotLimits();                       % get the limits for the plots we will make
  preStimMS = [trialStructs(:).preStimMS];                  	% get preStim times for each trial
  RTs = [trials(:).reactTimeMS];                              % get all trial RTs
  indices.early = (eotCodes == 1) & (preStimMS + RTs + plotRTStartMS > 200);
end
