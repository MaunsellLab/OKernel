function bootstraps = getCaseBootstraps(U, dataDirName, dataName, limits, mustRecreate)
  % Return an array containing the optogentic stimuli from the sessions in a table
  if nargin < 5
    mustRecreate = false;
  end
  if isempty(limits.oneDay)
    theName = strcat(dataDirName, ' Analysis/Mat Files/', dataName, {' '}, limits.animal, ' Profiles.mat');
  else
    theName = strcat(dataDirName, ' Analysis/Mat Files/', dataName, {' '}, limits.animal, {' '}, limits.oneDay, ' Profiles.mat');
  end
  bootFileName = theName{1};
  if isfile(bootFileName) && ~mustRecreate
    load(bootFileName, 'bootstraps');
  else
    bootstraps = struct('hitProfiles', [], 'missProfiles', [], 'RTProfiles', [], 'stimRTProfiles', [], ...
      'FAProfiles', []);
    for s = 1:size(U, 1)
      clear file trials dParams
      fileName = strcat(dataDirName, U.animal(s), '/MatFiles/', U.date(s));
      fprintf('%s\n', fileName);
      load(fileName);   %#ok<LOAD>
      if exist('trials', 'var') && isfield(trials, 'trial')
          row = U(s, :);
          row.rampMS = trials(1).trial.visualRampDurMS;
          [row.RTMinMS, row.RTMaxMS, row.missMinMS] = getRTParams(row.rampMS);
          bootstraps = getFileProfiles(bootstraps, file, trials, row);
      else
        fprintf('Error loading %s', fileName);
      end   
    end
    save(bootFileName, 'bootstraps');
  end
end

%%
function [row, hitIndices, missIndices, faIndices] = countsAndIndices(doStim, row, file, trials, eotCodes, stimIndices)

  trialStructs = [trials(:).trial];                           % move trialStructs to an array for access
	preStimMS = [trialStructs(:).preStimMS];                  	% get preStim times for each trial
  RTs = [trials(:).reactTimeMS];                              % get all trial RTs
  indices = stimIndices == doStim;                            % get indices for stim (or noStim) trials
	hitIndices = indices & eotCodes == 0;                       % get hit trials from OK
  faIndices = indices & eotCodes == 1;                        % get fa trials from OK
  missIndices = indices & eotCodes == 2;                      % get miss trials from OK
  
  % Calculate Hit & FA rates based simply on the EOT codes
  if doStim
    row.rawHits = sum(hitIndices) / (sum(hitIndices) + sum(missIndices)); % simple hit rate
    FAsPerMS = sum(faIndices) / sum(preStimMS(faIndices) + RTs(faIndices));
    row.rawFAs = 1.0 - (1.0 - FAsPerMS) ^ (file.rewardedLimitMS - file.tooFastMS);
  else
    row.rawNoStimHits = sum(hitIndices) / (sum(hitIndices) + sum(missIndices)); % simple hit rate
    FAsPerMS = sum(faIndices) / sum(preStimMS(faIndices) + RTs(faIndices));
    row.rawNoStimFAs = 1.0 - (1.0 - FAsPerMS) ^ (file.rewardedLimitMS - file.tooFastMS);
  end
  
  [~, ~, plotRTStartMS] = plotLimits();                       % get the limits for the plots we will make
  % for FAs, we must eliminate those that don't have enought time before the FA to avoid sampling before 
  % the end of the ramp (which is fixed at 200 ms). 
  faIndices = (indices & eotCodes == 1) & (preStimMS + RTs + plotRTStartMS > 200);
  earlyIndices = hitIndices & RTs < row.RTMinMS;              % hits outside of the allowed RT range
  lateIndices = hitIndices & RTs >= row.RTMaxMS;              % beyond latest RT counted as a hit
  faIndices = faIndices | earlyIndices;                       % add early hits into FAs
  missRTIndices = hitIndices & RTs >= row.missMinMS;          % hits so late they can add to misses
  hitIndices = hitIndices & ~(earlyIndices | lateIndices);    % strip out hits outside of range
  missIndices = (indices & eotCodes == 2) | missRTIndices;    % include very late hits into misses
  if doStim
    row.hits = sum(hitIndices);
    row.misses = sum(missIndices);
    row.FAs = sum(faIndices);
    row.numStim = row.hits + row.misses + row.numStim;        % don't count kEOTIgnored
  else
    row.noStimHits = sum(hitIndices);
    row.noStimMisses = sum(missIndices);
    row.noStimFAs = sum(faIndices); 
    row.numNoStim = row.noStimHits + row.noStimMisses + row.noStimFAs;  % don't count kEOTIgnored
  end
end

%%
function bootstraps = getFileProfiles(bootstraps, file, trials, row)

% Compute kernels for a session
  trialStructs = [trials(:).trial];
  meanPower = [trials(:).meanPowerMW];                        % get power applied for each trial                        
  if sum(meanPower) == 0                                      % no opto stimulation in this session
      return;
  end
  if isfield(trialStructs, 'pulseContrast')                   % get rid of any trials with reduced opto power
      stimIndices = meanPower > 0 & [trialStructs.pulseContrast] == 1;	% only trials with contrast == 1
  else
      stimIndices = meanPower > 0;
  end
  if sum(stimIndices) == 0                                    % no trials with full opto contrast
    return;
  end
  firstStimIndex = find(stimIndices > 0, 1);                 	% first stimulated trial
	lastStimIndex = find(stimIndices > 0, 1, 'last');          	% last stimulated trial
  
  % we are only going to consider the range between the first and last stimulated trials
  trials = trials(firstStimIndex:lastStimIndex);
%   trialStructs = trialStructs(firstStimIndex:lastStimIndex);
%   meanPower = meanPower(firstStimIndex:lastStimIndex);
  stimIndices = stimIndices(firstStimIndex:lastStimIndex);
  
  eotCodes = zeros(1, length(trials));                        % one, only one RT and endTrial per trial
  for t = 1:length(trials)                               
      if length(trials(t).reactTimeMS) > 1
          trials(t).reactTimeMS = trials(t).reactTimeMS(1);   
      end
      if ~isfield(trials(t), 'trialEnd') || isempty(trials(t).trialEnd)
          trials(t).trialEnd = -1;
      end
      eotCodes(t) = trials(t).trialEnd;
  end
  
  % get the counts for non-stimulated trials and then the stimulated trials
	[row, ~, ~, ~] = countsAndIndices(false, row, file, trials, eotCodes, stimIndices);
	[row, hitIndices, missIndices, faIndices] = countsAndIndices(true, row, file, trials, eotCodes, stimIndices);
  [plotStartMS, plotEndMS, plotRTStartMS] = plotLimits();     % get the limits for the plots we will make

  % get the hit kernel    
  hitProfiles = getStimProfiles(trials(hitIndices), plotStartMS, plotEndMS, true, false);
  bootstraps.hitProfiles = [bootstraps.hitProfiles; normProfiles(hitProfiles)];

  % get the RT aligned kernel
  RTProfiles = getStimProfiles(trials(hitIndices), plotRTStartMS, plotRTStartMS + plotEndMS - plotStartMS, true, true);
  bootstraps.RTProfiles = [bootstraps.RTProfiles; normProfiles(RTProfiles)];

  % get the Stim & RT aligned hit kernel
  stimRTProfiles = getStimRTProfiles(trials(hitIndices), plotStartMS, plotEndMS);
  bootstraps.stimRTProfiles = [bootstraps.stimRTProfiles; normProfiles(stimRTProfiles)];

  % get the miss kernel
  missProfiles = getStimProfiles(trials(missIndices), plotStartMS, plotEndMS, true, false);
  bootstraps.missProfiles = [bootstraps.missProfiles; normProfiles(missProfiles)];

  % Get the FA kernel
  if row.FAs > 0
      FAProfiles = getStimProfiles(trials(faIndices), plotRTStartMS, plotRTStartMS + plotEndMS - plotStartMS, true, true);
      bootstraps.FAProfiles = [bootstraps.FAProfiles; normProfiles(FAProfiles)];
  end
end

%%
function [RTMinMS, RTMaxMS, missMinMS] = getRTParams(rampMS)

  switch rampMS
    case 0
      RTMinMS = 200;
      RTMaxMS = 500;
      missMinMS = 750;
    case 500
      RTMinMS = 200;
      RTMaxMS = 500;
      missMinMS = 750;
    otherwise
      RTMinMS = 0;
      RTMaxMS = 0;
      missMinMS = 0;
  end
end

%%
function profiles = normProfiles(profiles)

meanPower = (max(max(profiles)) + min(min(profiles))) / 2.0;
profiles(profiles < meanPower) = -1.0;
profiles(profiles >= meanPower) = 1.0;
end
