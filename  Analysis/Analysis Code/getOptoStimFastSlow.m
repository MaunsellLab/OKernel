function [fileProfilesFast, fileProfilesSlow] = getOptoStimFastSlow(U, dataDirName)
  % Return an array containing the optogentic stimuli from the sessions in a table
  % Two arrays are returned: those from hit trials with RTs faster and slower than the median RT
  
  fileProfilesFast = [];
  fileProfilesSlow = [];
  for s = 1:size(U, 1)
    clear file trials dParams
    fileName = strcat(dataDirName, U.animal(s), '/MatFiles/', U.date(s));
    fprintf('%s\n', fileName);
    load(fileName, 'file', 'trials');
    if exist('trials', 'var') && isfield(trials, 'trial')
      row = U(s, :);
      row.rampMS = trials(1).trial.visualRampDurMS;
      [row.RTMinMS, row.RTMaxMS, row.missdMinMS] = getRTParams(row.rampMS);
      [fastProfiles, slowProfiles] = getFileProfiles(file, trials, row);
      fileProfilesFast = [fileProfilesFast; fastProfiles]; %#ok<*AGROW>
      fileProfilesSlow = [fileProfilesSlow; slowProfiles];
    else
      fprintf('Error loading %s', fileName);
    end   
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
function [fileProfilesFast, fileProfilesSlow] = getFileProfiles(file, trials, row)

% Return the opto stimuli from fast and slow trials separately
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
  
  % we are only going to consider the range between the first and last stimulated trials -- trim the vectors
  trials = trials(firstStimIndex:lastStimIndex);
  stimIndices = stimIndices(firstStimIndex:lastStimIndex);
%   trialStructs = trialStructs(firstStimIndex:lastStimIndex);
%   meanPower = meanPower(firstStimIndex:lastStimIndex);
  
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
  
  % get the counts for non-stimulated trials
	[row, ~, ~, ~] = countsAndIndices(false, row, file, trials, eotCodes, stimIndices);
  % get the counts for stimulated trials
	[~, hitIndices, ~, ~] = countsAndIndices(true, row, file, trials, eotCodes, stimIndices);
  [plotStartMS, plotEndMS, ~] = plotLimits();     % get the limits for the plots we will make

  % get the hit kernel.  There are no fast or slow misses
  allRTs = [trials(:).reactTimeMS];
  medianRT = median([trials(hitIndices).reactTimeMS]);
  fastHitIndices = hitIndices & allRTs <= medianRT;
  slowHitIndices = hitIndices & allRTs > medianRT;
  profiles = getStimProfiles(trials(fastHitIndices), plotStartMS, plotEndMS, true, false);
  fileProfilesFast = normProfiles(profiles);
  profiles = getStimProfiles(trials(slowHitIndices), plotStartMS, plotEndMS, true, false);
  fileProfilesSlow = normProfiles(profiles);
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
