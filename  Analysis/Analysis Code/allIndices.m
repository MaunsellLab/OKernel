function indices = allIndices(trials, eotCodes)
%
% Return logical arrays specifying valid indices for hit, fail and early trials.  The lists are modified
% according to trialType, which can select all, stimulated, or unstimulated trials.

%   selectIndices = true(1, length(stimIndices));               % take all trials
% 	indices.correct = selectIndices & eotCodes == 0;            % get hit trials
%   indices.early = selectIndices & eotCodes == 1;              % get early trials
%   indices.fail = selectIndices & eotCodes == 2;               % get fail trials
	indices.correct = eotCodes == 0;                            % get hit trials
  indices.early = eotCodes == 1;                              % get early trials
  indices.fail = eotCodes == 2;                               % get fail trials
  
  % for earlies, we must eliminate those that don't have enough time before the FA to avoid sampling before 
  % the end of the ramp (which is fixed at 200 ms). 
  trialStructs = [trials(:).trial];                           % extract trialStructs to an array for access
  preStimMS = [trialStructs(:).preStimMS];                  	% get preStim times for each trial
  
  % for the RT-aligned plot, we need to make sure that the plotRTStartMS before the RT doesn't 
  % go back into the 200 ms period where the optogenetic stimulus ramps up from zero.  If it does,
  % the kernels will all show the ramp
  [~, ~, plotRTStartMS] = plotLimits();                       % get the limits for the plots we will make
  RTs = [trials(:).reactTimeMS];                              % get all trial RTs
%   indices.early = (selectIndices & eotCodes == 1) & (preStimMS + RTs + plotRTStartMS > 200);
  indices.early = (eotCodes == 1) & (preStimMS + RTs + plotRTStartMS > 200);
end
