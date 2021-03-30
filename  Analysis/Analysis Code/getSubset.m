function [U, maxDate] = getSubset(mode, dataDirName, tableDataName, oneDay, limits)
  % Return a subset of files based on selection criteria.
  % If necessary, a new table of summary results will be compiled.
    
  [T, maxDate] = getTable(dataDirName, tableDataName);                % get table of all sessions
  switch mode
    case 'normal'
      controls = controlSessions(T);                                  % logic array flagging control sesions
      valid = T.rampMS == limits.rampMS & T.kernelCI > 0 & ~controls; % empty entries have zero for kernelCI
    case 'control'
       controls = controlTestSessions(T);                         	  % logic array flagging control test sesions
       valid = T.rampMS == limits.rampMS & T.kernelCI > 0 & controls; % empty entries have zero for kernelCI
   case 'prePostControl'
       controls = prePostControlSessions(T);                         	% logic array flagging control test sesions
       valid = T.rampMS == limits.rampMS & T.kernelCI > 0 & controls; % empty entries have zero for kernelCI
    otherwise
      fprintf('getSubset: unrecognized mode');
  end
  if ~strcmp(limits.animal, 'All')
    valid = valid & T.animal == limits.animal;
  end
  if ~isempty(oneDay)
    valid = valid & T.date == oneDay;
  end
  valid = valid & T.corrects > limits.minTrials;                        	% enough hit trials
  valid = valid & T.fails > limits.minTrials;                      	% enough fail trials
  valid = valid & T.kernelPeak ./ T.kernelCI > limits.criterion;
  if limits.minDec ~= -1
    stimHitRate = double(T.corrects) ./ double(T.corrects + T.fails);
    noStimHitRate = double(T.noStimCorrects) ./ double(T.noStimCorrects + T.noStimMisses);
    valid = valid & stimHitRate <= noStimHitRate - limits.minDec;     % min behavioral decrement in hit rate
  end
  U = T(valid, :);
  U = withoutTrainingDays(U);                                         % remove sessions from training periods
  
  % remove any animals with too few sessions
  if isempty(limits.oneDay) && limits.minSessions > 0
    animals = unique(U.animal);
    for a = 1:length(animals)
      animalRows = U.animal == animals{a};
      if sum(animalRows) < limits.minSessions
        U = U(~animalRows, :);
      end
    end
  end
  
  if size(U, 1) == 0
    fprintf('getSubset: No valid sessions found for %s rampMS %d, animal %s minTrials %d and decrement %.2f\n', ...
      mode, limits.rampMS, limits.animal, limits.minTrials, limits.minDec);
  end
end

%%
function controls = controlSessions(T)

% control measurement sesssion.  These are contralateral stimulation
% controls that should be excluded from general analysis

  cSessions = {
    {'1218', {'2020-06-06', '2020-06-07', '2020-06-08', '2020-06-09', '2020-06-10', '2020-06-11', '2020-06-12', ...
              '2020-06-13', '2020-06-14', '2020-06-15', '2020-06-16'}},...
    {'1220', {'2020-06-22', '2020-06-23', '2020-06-24', '2020-06-25', '2020-06-26', '2020-06-27', '2020-06-28'}},...
    {'1257', {'2020-05-30', '2020-05-31', '2020-06-01', '2020-06-02', '2020-06-03', '2020-06-04', '2020-06-05', ...
              '2020-06-06', '2020-06-07', '2020-06-08', '2020-06-09', '2020-06-10', '2020-06-11', '2020-06-12', ...
              '2020-06-13', '2020-06-14'}},...
  };
  
  controls = false(height(T), 1);
  for a = 1:length(cSessions)
    for d = 1:length(cSessions{a}{2})
      controls = controls | (T.animal == cSessions{a}{1} & T.date == cSessions{a}{2}{d});
    end
  end
end

function controls = prePostControlSessions(T)

% control measurement sesssion.  These are a subset of normal stimulation
% sessions that were done immediately before and after the control stimulation
% session.  They are matched in number to the control stimulation sessions
% to keep the S/N balanced between the two

  cSessions = {
    {'1218', {'2020-06-01', '2020-06-02', '2020-06-03', '2020-06-04', '2020-06-05', ...
              '2020-06-21', '2020-06-22', '2020-06-23', '2020-06-24', '2020-06-25'}},...
    {'1220', {'2020-06-17', '2020-06-18', '2020-06-19', '2020-06-20', '2020-06-21'}},...
    {'1257', {'2020-05-25', '2020-05-26', '2020-05-27', '2020-05-28', '2020-05-29', ...
              '2020-06-21', '2020-06-22', '2020-06-23', '2020-06-24', '2020-06-25'}},...
  };
  
  controls = false(height(T), 1);
  for a = 1:length(cSessions)
    for d = 1:length(cSessions{a}{2})
      controls = controls | (T.animal == cSessions{a}{1} & T.date == cSessions{a}{2}{d});
    end
  end
end

function controls = controlTestSessions(T)

% control measurement sesssion. These are contralateral stimulation
% controls that should be excluded from general analysis.  This list is
% balanced to have the same number of sessions as the pre/post control
% measurements, to keep the S/N balanced between the two.

  cSessions = {
    {'1218', {'2020-06-07', '2020-06-08', '2020-06-09', '2020-06-10', '2020-06-11', ...
              '2020-06-12', '2020-06-13', '2020-06-14', '2020-06-15', '2020-06-16'}},...
    {'1220', {'2020-06-22', '2020-06-23', '2020-06-24', '2020-06-25', '2020-06-28'}},...
    {'1257', {'2020-06-04', '2020-06-05', '2020-06-06', '2020-06-07', '2020-06-09', ...
              '2020-06-10', '2020-06-11', '2020-06-12', '2020-06-13', '2020-06-14'}},...
  };
  
  controls = false(height(T), 1);
  for a = 1:length(cSessions)
    for d = 1:length(cSessions{a}{2})
      controls = controls | (T.animal == cSessions{a}{1} & T.date == cSessions{a}{2}{d});
    end
  end
end

%%
% function [row, indices] = countsAndIndices(doStim, row, file, trials, eotCodes, stimIndices)
% 
%   trialStructs = [trials(:).trial];                           % move trialStructs to an array for access
% 	preStimMS = [trialStructs(:).preStimMS];                  	% get preStim times for each trial
%   RTs = [trials(:).reactTimeMS];                              % get all trial RTs
%   allIndices = stimIndices == doStim;                         % get indices for stim (or noStim) trials
% 	indices.correct = allIndices & eotCodes == 0;               % get hit trials
%   indices.early = allIndices & eotCodes == 1;                    % get fa trials
%   indices.fail = allIndices & eotCodes == 2;                  % get fail trials
%   
%   % Calculate Hit & FA rates based simply on the EOT codes
%   if doStim
%     row.rawCorrects = sum(indices.correct) / (sum(indices.correct) + sum(indices.fail)); % simple hit rate
%     earliesPerMS = sum(indices.early) / sum(preStimMS(indices.early) + RTs(indices.early));
%     row.rawEarlies = 1.0 - (1.0 - earliesPerMS) ^ (file.rewardedLimitMS - file.tooFastMS);
%   else
%     row.rawNoStimCorrects = sum(indices.correct) / (sum(indices.correct) + sum(indices.fail)); % simple hit rate
%     earliesPerMS = sum(indices.early) / sum(preStimMS(indices.early) + RTs(indices.early));
%     row.rawNoStimEarlies = 1.0 - (1.0 - earliesPerMS) ^ (file.rewardedLimitMS - file.tooFastMS);
%   end
%   
%   [~, ~, plotRTStartMS] = plotLimits();                       % get the limits for the plots we will make
%   % for earlies, we must eliminate those that don't have enough time before the early to avoid sampling before 
%   % the end of the ramp (which is fixed at 200 ms). 
%   indices.early = (allIndices & eotCodes == 1) & (preStimMS + RTs + plotRTStartMS > 200);
%   earlyIndices = indices.correct & RTs < row.RTMinMS;             % corrects outside of the allowed RT range
%   lateIndices = indices.correct & RTs >= row.RTMaxMS;             % beyond latest RT counted as a hit
%   indices.early = indices.early | earlyIndices;                   % add early corrects into earlies
%   failRTIndices = indices.correct & RTs >= row.failMinMS;         % corrects so late they can add to failes
%   indices.correct = indices.correct & ~(earlyIndices | lateIndices);  % strip out corrects outside of range
%   indices.fail = (allIndices & eotCodes == 2) | failRTIndices;    % include very late corrects into failes
%   if doStim
%     row.corrects = sum(indices.correct);
%     row.fails = sum(indices.fail);
%     row.earlies = sum(indices.early);
%     row.numStim = row.corrects + row.fails + row.numStim;        % don't count kEOTIgnored
%   else
%     row.noStimCorrects = sum(indices.correct);
%     row.noStimMisses = sum(indices.fail);
%     row.noStimEarlies = sum(indices.early); 
%     row.numNoStim = row.noStimCorrects + row.noStimMisses + row.noStimEarlies;  % don't count kEOTIgnored
%   e% function [row, indices] = countsAndIndices(doStim, row, file, trials, eotCodes, stimIndices)
% 
%   trialStructs = [trials(:).trial];                           % move trialStructs to an array for access
% 	preStimMS = [trialStructs(:).preStimMS];                  	% get preStim times for each trial
%   RTs = [trials(:).reactTimeMS];                              % get all trial RTs
%   allIndices = stimIndices == doStim;                         % get indices for stim (or noStim) trials
% 	indices.correct = allIndices & eotCodes == 0;               % get hit trials
%   indices.early = allIndices & eotCodes == 1;                 % get early trials
%   indices.fail = allIndices & eotCodes == 2;                  % get fail trials
%   
%   % Calculate Hit & FA rates based simply on the EOT codes
%   if doStim
%     row.rawCorrects = sum(indices.correct) / (sum(indices.correct) + sum(indices.fail)); % simple hit rate
%     earliesPerMS = sum(indices.early) / sum(preStimMS(indices.early) + RTs(indices.early));
%     row.rawEarlies = 1.0 - (1.0 - earliesPerMS) ^ (file.rewardedLimitMS - file.tooFastMS);
%   else
%     row.rawNoStimCorrects = sum(indices.correct) / (sum(indices.correct) + sum(indices.fail)); % simple hit rate
%     earliesPerMS = sum(indices.early) / sum(preStimMS(indices.early) + RTs(indices.early));
%     row.rawNoStimEarlies = 1.0 - (1.0 - earliesPerMS) ^ (file.rewardedLimitMS - file.tooFastMS);
%   end
%   
%   [~, ~, plotRTStartMS] = plotLimits();                       % get the limits for the plots we will make
%   % for earlies, we must eliminate those that don't have enough time before the FA to avoid sampling before 
%   % the end of the ramp (which is fixed at 200 ms). 
%   indices.early = (allIndices & eotCodes == 1) & (preStimMS + RTs + plotRTStartMS > 200);
%   earlyIndices = indices.correct & RTs < row.RTMinMS;             % corrects outside of the allowed RT range
%   lateIndices = indices.correct & RTs >= row.RTMaxMS;             % beyond latest RT counted as a hit
%   indices.early = indices.early | earlyIndices;                         % add early corrects into earlies
%   failRTIndices = indices.correct & RTs >= row.failMinMS;         % corrects so late they can add to fails
%   indices.correct = indices.correct & ~(earlyIndices | lateIndices);  % strip out corrects outside of range
%   indices.fail = (allIndices & eotCodes == 2) | failRTIndices;    % include very late corrects into fails
%   if doStim
%     row.corrects = sum(indices.correct);
%     row.fails = sum(indices.fail);
%     row.earlies = sum(indices.early);
%     row.numStim = row.corrects + row.fails + row.numStim;        % don't count kEOTIgnored
%   else
%     row.noStimCorrects = sum(indices.correct);
%     row.noStimMisses = sum(indices.fail);
%     row.noStimEarlies = sum(indices.early); 
%     row.numNoStim = row.noStimCorrects + row.noStimMisses + row.noStimEarlies;  % don't count kEOTIgnored
%   end
% end
% end

%%
function fileNames = getFileNames(dirName)

    dirStructs = dir(dirName);                              % data directory contents
    fileNames = {dirStructs(:).name};                       % data directory file names
    numFiles = length(fileNames);
    validFiles = false(1, numFiles);                        % find only files with number names (animals)
    for f = 1:numFiles
        validFiles(f) = isstrprop(fileNames{f}(1), 'digit') && length(fileNames{f}) > 1;
    end
    fileNames = {fileNames{validFiles}}; %#ok<*CCAT1>
    fileValues = str2double(fileNames);             % sort the files numerically
    [~, indices] = sort(fileValues);
    fileNames = {fileNames{indices}};
end

function indices = getIndices(row, trialType, earlyType, trials, eotCodes, stimIndices)
%
% Return logical arrays specifying valid indices for hit, fail and early trials.  The lists are modified
% according to trialType, which can select all, stimulated, or unstimulated trials.

  trialStructs = [trials(:).trial];                           % extract trialStructs to an array for access
  RTs = [trials(:).reactTimeMS];                              % get all trial RTs
  trialType = lower(trialType);
  if strcmp(trialType, 'stim') || strcmp(trialType, 'nostim')
    selectIndices = stimIndices == strcmp(trialType, 'stim');	% get logical array for for stim (or noStim) trials
  else
     selectIndices = true(1, length(stimIndices));            % take all trials
  end
	indices.correct = selectIndices & eotCodes == 0;            % get hit trials
  indices.early = selectIndices & eotCodes == 1;              % get early trials
  indices.fail = selectIndices & eotCodes == 2;               % get fail trials
  
  % for earlies, we must eliminate those that don't have enough time before the FA to avoid sampling before 
  % the end of the ramp (which is fixed at 200 ms). 
  if strcmp(earlyType, 'corrected')
    preStimMS = [trialStructs(:).preStimMS];                  	% get preStim times for each trial
    [~, ~, plotRTStartMS] = plotLimits();                       % get the limits for the plots we will make
    indices.early = (selectIndices & eotCodes == 1) & (preStimMS + RTs + plotRTStartMS > 200);
    earlyIndices = indices.correct & RTs < row.RTMinMS;        	% corrects outside of the allowed RT range
    lateIndices = indices.correct & RTs >= row.RTMaxMS;         % beyond latest RT counted as a hit
    failRTIndices = indices.correct & RTs >= row.failMinMS;    	% corrects so late they can add to fails
    indices.early = indices.early | earlyIndices;                    	% add early corrects into earlies
    indices.correct = indices.correct & ~(earlyIndices | lateIndices);  % strip out corrects outside of range
    indices.fail = (selectIndices & eotCodes == 2) | failRTIndices;    % include very late corrects into fails
  end
end

%%
function row = getKernels(file, trials, row)

% Compute kernels for a session
  trialStructs = [trials(:).trial];
  meanPower = [trials(:).meanPowerMW];                        % get power applied for each trial                        
  if sum(meanPower) < 10                                      % no opto stimulation in this session
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
  trialStructs = trialStructs(firstStimIndex:lastStimIndex);
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
  
	preStimMS = [trialStructs(:).preStimMS];                  	% get preStim times for each trial  
  RTs = [trials(:).reactTimeMS];                              % get all trial RTs
  % Add to the row the counts for non-stimulated trials
% 	[row, ~] = countsAndIndices(false, row, file, trials, eotCodes, stimIndices);
  indices = getIndices(row, 'nostim', 'uncorrected', trials, eotCodes, stimIndices);
  row.rawNoStimCorrects = sum(indices.correct) / (sum(indices.correct) + sum(indices.fail)); % simple hit rate
  earliesPerMS = sum(indices.early) / sum(preStimMS(indices.early) + RTs(indices.early));
  row.rawNoStimEarlies = 1.0 - (1.0 - earliesPerMS) ^ (file.rewardedLimitMS - file.tooFastMS);
  indices = getIndices(row, 'nostim', 'corrected', trials, eotCodes, stimIndices);
  row.noStimCorrects = sum(indices.correct);
  row.noStimMisses = sum(indices.fail);
  row.noStimEarlies = sum(indices.early); 
  row.numNoStim = row.noStimCorrects + row.noStimMisses + row.noStimEarlies;  % don't count kEOTIgnored

  % add to the row counts for stimulated trials
% 	[row, indices] = countsAndIndices(true, row, file, trials, eotCodes, stimIndices);
  indices = getIndices(row, 'stim', 'uncorrected', trials, eotCodes, stimIndices);
  row.rawCorrects = sum(indices.correct) / (sum(indices.correct) + sum(indices.fail)); % simple hit rate
  earliesPerMS = sum(indices.early) / sum(preStimMS(indices.early) + RTs(indices.early));
  row.rawEarlies = 1.0 - (1.0 - earliesPerMS) ^ (file.rewardedLimitMS - file.tooFastMS);
  indices = getIndices(row, 'stim', 'corrected', trials, eotCodes, stimIndices);
  row.corrects = sum(indices.correct);
  row.fails = sum(indices.fail);
  row.earlies = sum(indices.early);
  row.numStim = row.corrects + row.fails + row.numStim;                  % don't count kEOTIgnored

  % calculate the overall d'
  allIndices = getIndices(row, 'all', 'corrected', trials, eotCodes, stimIndices);
  % If we don't have indices at this point, there's nothing left to do.
  if sum(allIndices.correct | allIndices.fail | allIndices.early) == 0
    return;
  end
	[respLimitsMS, allIndices, ~, ~] = getResponseLimits(file, trials, allIndices);
  hitRate = sum(allIndices.correct) / (sum(allIndices.correct) + sum(allIndices.fail));
  earlyRate = earlyRate(file, trials, allIndices.correct, allIndices.fail, allIndices.early);
  row.RTWindowMS = diff(respLimitsMS);
  row.pFA = 1.0 - exp(-earlyRate * row.RTWindowMS / 1000.0);
  row.pHit = (hitRate - row.pFA) / (1.0 - row.pFA);
  [row.dPrime, row.c] = dprime(row.pHit, row.pFA);
%     fprintf('  All:     pFA = %.2f, pHit = %.2f d-prime = %.1f, c = %.1f\n', row.pFA, row.pHit, row.dPrime, row.c);
  
  % calculate the nostim trial d'
%   indices = getIndices(row, 'nostim', 'corrected', trials, eotCodes, stimIndices);
% 	[respLimitsMS, indices, ~, ~] = getResponseLimits(file, trials, indices);
  indices.correct = allIndices.correct & ~stimIndices;
  indices.fail = allIndices.fail & ~stimIndices;
  indices.early = allIndices.early & ~stimIndices;
  if sum(indices.correct | indices.fail | indices.early) > 0
    hitRate = sum(indices.correct) / (sum(indices.correct) + sum(indices.fail));
    earlyRate = earlyRate(file, trials, indices.correct, indices.fail, indices.early);
    row.noStimPFA = 1.0 - exp(-earlyRate * row.RTWindowMS / 1000.0);
    row.noStimPHit = (hitRate - row.pFA) / (1.0 - row.pFA);         % using overall pFA
    [row.noStimDPrime, row.noStimC] = dprime(row.noStimPHit, row.pFA);
    fprintf('  No Stim: pFA = %.2f, pHit = %.2f d-prime = %.1f, c = %.1f\n', ...
              row.noStimPFA, row.noStimPHit, row.noStimDPrime, row.noStimC);
  end
    
  % calculate the stim trial d'
%   indices = getIndices(row, 'stim', 'corrected', trials, eotCodes, stimIndices);
% 	[respLimitsMS, indices, ~, ~] = getResponseLimits(file, trials, indices);
  indices.correct = allIndices.correct & stimIndices;
  indices.fail = allIndices.fail & stimIndices;
  indices.early = allIndices.fa & stimIndices;
  if sum(indices.correct | indices.fail | indices.early) > 0
    hitRate = sum(indices.correct) / (sum(indices.correct) + sum(indices.fail));
    earlyRate = earlyRate(file, trials, indices.correct, indices.fail, indices.early);
    row.stimPFA = 1.0 - exp(-earlyRate * row.RTWindowMS / 1000.0);
    row.stimPHit = (hitRate - row.pFA) / (1.0 - row.pFA);   % using overall pFA
    [row.stimDPrime, row.stimC] = dprime(row.stimPHit, row.pFA);
    fprintf('  Stim:    pFA = %.2f, pHit = %.2f d-prime = %.1f, c = %.1f  delta: %.5f\n', ...
              row.stimPFA, row.stimPHit, row.stimDPrime, row.stimC, row.noStimDPrime - row.stimDPrime);
  end
  
  % If we don't have indices at this point, there's nothing left to do.
  
  if sum(indices.correct | indices.fail | indices.early) == 0
    return;
  end
  
  [plotStartMS, plotEndMS, plotRTStartMS] = plotLimits();     % get the limits for the plots we will display  
  % get the hit kernel.  We use the indices from the last call to getResponseLimits, which will be all the 
  % stimulated and unstimulated trials within the detected response window.
  profiles = getStimProfiles(trials(indices.correct), plotStartMS, plotEndMS, true, false);
  hitSums(:) = normSum(profiles);
  hitKernel = hitSums / row.corrects;
  row.hitKernel = {hitKernel};
  hitCI = stimCI(row.corrects);

  % get the RT aligned kernel
  profiles = getStimProfiles(trials(indices.correct), plotRTStartMS, plotRTStartMS + plotEndMS - plotStartMS, true, true);
  RTSum(:) = normSum(profiles);
  row.RTKernel = {RTSum / row.corrects};

  % get the Stim & RT aligned hit kernel
  profiles = getStimRTProfiles(trials(indices.correct), plotStartMS, plotEndMS);
  SRTSum(:) = normSum(profiles);
  row.SRTKernel = {SRTSum / row.corrects};

  % get the fail kernel
  profiles = getStimProfiles(trials(indices.fail), plotStartMS, plotEndMS, true, false);
  failSum(:) = normSum(profiles);
  failKernel = failSum / row.fails;
  row.failKernel = {failKernel};
  failCI = stimCI(row.fails);

  % Get the FA kernel
  if row.earlies > 0
      profiles = getStimProfiles(trials(indices.early), plotRTStartMS, plotRTStartMS + plotEndMS - plotStartMS, true, true);
      FASum(:) = normSum(profiles);
      row.earlies = size(profiles, 1);                  % getStimProfiles might reject some trials as too short
      row.FAKernel = {FASum / row.earlies};
  end

  % If we have corrects and fails, save the hit, fail, total and random kernels.
  row.kernelCI = sqrt(hitCI^2 + failCI^2);
  if ~isempty(hitKernel) && ~isempty(failKernel)
    row.peakMinMS = plotStartMS;
    row.peakMaxMS = plotEndMS;
    startIndex = row.peakMinMS - plotStartMS + 1;
    endIndex = row.peakMaxMS - plotStartMS;
    row.kernelPeak =  max(abs((hitKernel(startIndex:endIndex) - failKernel(startIndex:endIndex))) / row.kernelCI);
    row.randomKernel = getRandomKernel(row.corrects, row.fails, trialStructs(1).pulseDurMS, plotEndMS - plotStartMS);
  end

  % add to the RT distributions
  row.correctRTs = {[trials([trials(:).trialEnd] == 0).reactTimeMS]};
  row.wrongRTs = {[trials([trials(:).trialEnd] == 1).reactTimeMS]};
  failRTs = [trials([trials(:).trialEnd] == 2).reactTimeMS];
  failRTs(failRTs < 0) = 100000;        % include fails in count, but don't let them display on plot
  failRTs(failRTs > 100000) = 100000;
  row.failRTs = {failRTs};
end

%% getTable -- load or create the table with all the summary data
% function [T, maxDate] = getTable(dataDirName, tableDataName, peakMinMS, peakMaxMS, plotStartMS, plotEndMS,...
%   plotRTStartMS)
function [T, maxDate] = getTable(dataDirName, tableDataName)

  if ~isfile(tableDataName)
      varNames = {'animal', 'date', 'numStim', 'corrects', 'fails', 'earlies', ...
          'numNoStim', 'noStimCorrects', 'noStimMisses', 'noStimEarlies', ...
          'rawCorrects', 'rawEarlies', 'rawNoStimCorrects', 'rawNoStimEarlies' ...
          'RTWindowMS', 'pHit', 'pFA', 'dPrime', 'c', ...
          'stimPHit', 'stimPFA', 'stimDPrime', 'stimC', ...
          'noStimPHit', 'noStimPFA', 'noStimDPrime', 'noStimC', ...
          'RTMinMS', 'RTMaxMS', ...
          'failMinMS', 'peakMinMS', 'peakMaxMS', 'rampMS', 'kernelCI', 'kernelPeak', ...
          'hitKernel', 'failKernel', 'RTKernel', 'FAKernel', 'SRTKernel', 'randomKernel', 'correctRTs', 'failRTs', ...
          'wrongRTs'};
      varTypes = {'string', 'string', 'uint32', 'uint32', 'uint32', 'uint32', ...
          'uint32', 'uint32', 'uint32', 'uint32', ...
          'double', 'double', 'double', 'double', ...
          'double', 'double', 'double',  'double', 'double', ...
          'double', 'double', 'double',  'double', ...
          'double', 'double', 'double',  'double', ...
          'uint32', 'uint32', ...
          'uint32', 'double', 'double', 'uint32', 'double', 'double', ...
          'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', ...
          'cell'};
      T = table('size', [0, length(varNames)], 'variableTypes', varTypes, 'variableNames', varNames);
      save(tableDataName, 'T');
      fprintf('No data table found: Created one.\n');
  end
  load(tableDataName, 'T');
  % count the number of sessions
  animalNames = getFileNames(dataDirName);
  numAnimals = length(animalNames);
  numSessions = 0;
  for a = 1:numAnimals
      fileNames = getFileNames([dataDirName animalNames{a} '/MatFiles/']);
      numSessions = numSessions + length(fileNames);
  end
  % check every session, compile any missing entries calling getRowEntry().
  session = 1;
  maxDate = 0;
  for a = 1:numAnimals
      fileNames = getFileNames([dataDirName animalNames{a} '/MatFiles/']);
      numFiles = length(fileNames);
      for f = 1:numFiles                                  	% for each file
          [~, fileName, ~] = fileparts(fileNames{f});      	% get file name
          if length(fileName) > 10                        	% other sorts of data in folder
             continue; 
          end
          maxDate = max(maxDate, str2double(strrep(fileName(3:end), '-', '')));
          row = getTableRow(T, animalNames{a}, fileName);
          if isempty(row)
              rowEntry.animal = animalNames{a};
              rowEntry.date = fileName;
              rowEntry.peakMinMS = -400;
              rowEntry.peakMaxMS = 0;
              rowEntry = getRowEntry(dataDirName, rowEntry);
              if isempty(rowEntry)
                disp(row);
              end
              T = [T; struct2table(rowEntry, 'asArray', 1)];   %#ok<AGROW>
              fprintf('Session %d of %d: %s %s\n', session, numSessions, animalNames{a}, fileName);
          end
          session = session + 1;
      end
  end
  save(tableDataName, 'T'); 
end

%%
function row = getRowEntry(dataDirName, row)

    % initialize to an empty row
    row.numStim = 0;
    row.corrects = 0;
    row.fails = 0;
    row.earlies = 0;
    row.numNoStim = 0;
    row.noStimCorrects = 0;
    row.noStimMisses = 0;
    row.noStimEarlies = 0;
    row.rawCorrects = 0.0;
    row.rawEarlies = 0.0;
    row.rawNoStimCorrects = 0.0;
    row.rawNoStimEarlies = 0.0;
    row.pHit = 0.0;
    row.pFA = 0.0;
    row.dPrime = 0.0;
    row.c = 0.0;
    row.RTWindowMS = 0.0;
    row.stimPHit = 0.0;
    row.stimPFA = 0.0;
    row.stimDPrime = 0.0;
    row.stimC = 0.0;
    row.noStimPHit = 0.0;
    row.noStimPFA = 0.0;
    row.noStimDPrime = 0.0;
    row.noStimC = 0.0;
    row.kernelCI = 0.0;
    row.kernelPeak = 0.0;
    row.hitKernel = {0};
    row.failKernel = {0};
    row.RTKernel = {0};
    row.SRTKernel = {0};
    row.FAKernel = {0};
    row.randomKernel = {0};
    row.correctRTs = {0};
    row.failRTs = {0};
    row.wrongRTs = {0};
    row.rampMS = 0;                                           % no data for the requested row
    row.RTMinMS = 0;
    row.RTMaxMS = 0;
    row.failMinMS = 0;
    % load the data
    clear file trials dParams
    load([dataDirName row.animal '/MatFiles/' row.date]);     %#ok<LOAD>
    if exist('trials', 'var') && isfield(trials, 'trial')     %#ok<NODEF>
        row.rampMS = trials(1).trial.visualRampDurMS;
        [row.RTMinMS, row.RTMaxMS, row.failMinMS] = getRTParams(row.rampMS);
        trials = validateTrials(row, trials);
        row = getKernels(file, trials, row);    
    end
end

%%
function [RTMinMS, RTMaxMS, failMinMS] = getRTParams(rampMS)

  switch rampMS
    case 0
      RTMinMS = 200;
      RTMaxMS = 500;
      failMinMS = 750;
    case 500
      RTMinMS = 200;
      RTMaxMS = 500;
      failMinMS = 750;
    otherwise
      RTMinMS = 200;
      RTMaxMS = 500;
      failMinMS = 750;
  end
end

%%
function row = getTableRow(T, name, date)
    
    valid = strcmp(T.animal(:), name) & strcmp(T.date(:), date);
    row = find(valid);
    if length(row) > 1                                          % if there are multiple entries, delete all but first
        T(:,row(2:end)) = [];                                   %#ok<NASGU>
        row = row(1);
    end
end

%%
function normSums = normSum(profiles)

meanPower = (max(max(profiles)) + min(min(profiles))) / 2.0;
profiles(profiles < meanPower) = 0;
profiles(profiles >= meanPower) = 1;
normSums = sum(profiles, 1);                                    % summed normalized profiles
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
