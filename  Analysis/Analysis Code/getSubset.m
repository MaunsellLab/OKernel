function U = getSubset(mode, dataDirName, tableDataName, limits)
  % Return a subset of files based on selection criteria.
  % If necessary, a new table of summary results will be compiled.
    
  T = getTable(dataDirName, [dataDirName, tableDataName]);                % get table of all sessions
  switch mode
    case 'normal'
      controls = controlSessions(T);                                  % logic array flagging control sesions
      valid = any(T.rampMS == limits.rampMS & T.kernelCI > 0, 2) & ~controls; % empty entries have zero for kernelCI
    otherwise
      fprintf('getSubset: unrecognized mode');
  end
  if ~strcmp(limits.animal, 'All')
    valid = valid & T.animal == limits.animal;
  end
  if ~isempty(limits.oneDay)
    valid = valid & T.date == limits.oneDay;
  end
  valid = valid & (T.numStim > 0);                                    % only looking at sessions with stimulation
  valid = valid & (T.stimCorrects > limits.minTrials);                % enough correct trials
  valid = valid & T.stimFails > limits.minTrials;                    	% enough fail trials
  valid = valid & T.kernelPeak > limits.criterion;
  
  % d' no stim and decrement limits
  if limits.minDPrime ~= -1
    valid = valid & double(T.noStimDPrime) >= limits.minDPrime;
  end
  if limits.minDec ~= -1
    valid = valid & double(T.noStimDPrime) - double(T.stimDPrime) >= limits.minDec;     % min behavioral decrement in hit rate
  end
  if limits.maxMeanPowerMW ~= -1
    valid = valid & T.meanPowerMW <= limits.maxMeanPowerMW;
  end
  U = T(valid, :);
 	U.dPrime(U.dPrime == Inf) = NaN; 
	U.noStimDPrime(U.noStimDPrime == Inf) = NaN; 
	U.stimDPrime(U.stimDPrime == Inf) = NaN; 
 
  % Some limits related to over-session performance.  We can requie a minimum number of sessions for each animal/ramp,
  % and a minimum average delta-d'.
  if isempty(limits.oneDay)
    animals = unique(U.animal);
    rampDurs = unique(limits.rampMS);
    for r = 1:length(rampDurs)
      for a = 1:length(animals)
        animalRows = U.animal == animals{a} & U.rampMS == rampDurs(r);
        if limits.minSessions > 0 && sum(animalRows) < limits.minSessions
          U = U(~animalRows, :);
          continue;
        end
        if limits.minAvgDeltaDPrime >= 0 && nanmean(U.noStimDPrime(animalRows) - U.stimDPrime(animalRows)) < limits.minAvgDeltaDPrime
          fprintf('%s mean too low %.2f (n=%d)\n', ...
              animals{a}, nanmean(U.noStimDPrime(animalRows) - U.stimDPrime(animalRows)), sum(animalRows));
          U = U(~animalRows, :);
        else
          fprintf('%s mean is good %.2f (n=%d)\n', ...
              animals{a}, nanmean(U.noStimDPrime(animalRows) - U.stimDPrime(animalRows)), sum(animalRows));
        end
      end
    end
  end
  if size(U, 1) == 0
    if length(limits.rampMS) == 1
      fprintf('getSubset: No valid sessions found for ''%s'' for animal ''%s'' on %d ms ramps, minTrials %d and decrement %.2f\n', ...
        mode, limits.animal, limits.rampMS, limits.minTrials, limits.minDec);
    else
      fprintf('getSubset: No valid sessions found for ''%s'' for animal ''%s'' on multiple ramps, minTrials %d and decrement %.2f\n', ...
        mode, limits.animal, limits.minTrials, limits.minDec);
    end
  end
end

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

%%
function row = getKernels(file, trials, row)
%
% Compute kernels for one session
%
  trialStructs = [trials(:).trial];
  meanPower = [trials(:).meanPowerMW];                        % get power applied for each trial 
  row.meanPowerMW = mean(meanPower);
  row.maxPowerMW = max(meanPower);
  if isfield(trialStructs, 'pulseContrast')                   % get rid of any trials with reduced opto power
      stimIndices = meanPower > 0 & [trialStructs.pulseContrast] == 1;	% only trials with contrast == 1
  else
      stimIndices = meanPower > 0;
  end
  if sum(stimIndices) == 0                                    % no stimulation with full opto contrast
    return;
  end
  
  % we only consider the range between the first and last stimulated trials
  firstStimIndex = find(stimIndices > 0, 1);                 	% first stimulated trial
	lastStimIndex = find(stimIndices > 0, 1, 'last');          	% last stimulated trial
  trials = trials(firstStimIndex:lastStimIndex);
  trialStructs = trialStructs(firstStimIndex:lastStimIndex);
  stimIndices = stimIndices(firstStimIndex:lastStimIndex);
  
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
  
  [plotStartMS, plotEndMS, plotRTStartMS] = plotLimits();     % get the limits for the plots we will display  
  % get the hit kernel.  We use the indices set up to have only stimulated trials within the detected response window.
  profiles = getStimProfiles(trials(indices.correct), plotStartMS, plotEndMS, true, false);
%   hitSums(:) = normSum(profiles);
  hitSums(:) = sum(profiles, 1);
  hitKernel = hitSums / row.stimCorrects;
  row.hitKernel = {hitKernel};
  hitCI = stimCI(row.stimCorrects);

  % get the RT aligned kernel
  profiles = getStimProfiles(trials(indices.correct), plotRTStartMS, plotRTStartMS + plotEndMS - plotStartMS, true, true);
%   RTSum(:) = normSum(profiles);
  RTSum(:) = sum(profiles, 1);
  row.RTKernel = {RTSum / row.stimCorrects};

  % get the Stim & RT aligned hit kernel
  profiles = getStimRTProfiles(trials(indices.correct), plotStartMS, plotEndMS);
%   SRTSum(:) = normSum(profiles);
  SRTSum(:) = sum(profiles, 1);
  row.SRTKernel = {SRTSum / row.stimCorrects};

  % get the fail kernel
  profiles = getStimProfiles(trials(indices.fail), plotStartMS, plotEndMS, true, false);
  failSum(:) = sum(profiles, 1);
  failKernel = failSum / row.stimFails;
  row.failKernel = {failKernel};
  failCI = stimCI(row.stimFails);

  % Get the early kernel.  Eliminate trials that start before the end of the ramping stimulus.
	earlyIndices = indices.early & (preStimMS + RTs + plotRTStartMS > 200);
  if row.stimEarlies > 0
      profiles = getStimProfiles(trials(earlyIndices), plotRTStartMS, plotRTStartMS + plotEndMS - plotStartMS, true, true);
      earlySum(:) = normSum(profiles);
      row.stimEarlies = size(profiles, 1);                  % getStimProfiles might reject some trials as too short
      row.earlyKernel = {earlySum / row.stimEarlies};
  end

  % If we have corrects and fails, save the hit, fail, total and random kernels.
  row.kernelCI = sqrt(hitCI^2 + failCI^2);
  if ~isempty(hitKernel) && ~isempty(failKernel)
    row.peakMinMS = plotStartMS;
    row.peakMaxMS = plotEndMS;
    startIndex = row.peakMinMS - plotStartMS + 1;
    endIndex = row.peakMaxMS - plotStartMS;
    row.kernelPeak =  max(abs((hitKernel(startIndex:endIndex) - failKernel(startIndex:endIndex))) / row.kernelCI);
    row.randomKernel = getRandomKernel(row.stimCorrects, row.stimFails, trialStructs(1).pulseDurMS, plotEndMS - plotStartMS);
  end

  % add to the RT distributions
  row.correctRTs = {[trials(theIndices.correct).reactTimeMS]};
  row.earlyRTs = {[trials(theIndices.early).reactTimeMS]};
  failRTs = [trials(theIndices.fail).reactTimeMS];
  failRTs(failRTs < 0 | failRTs > 100000) = 100000;     % include fails in count, but don't let them display on plot
  row.failRTs = {failRTs};
end

%% getTable -- load or create the table with all the summary data
% function [T, maxDate] = getTable(dataDirName, tableDataName, peakMinMS, peakMaxMS, plotStartMS, plotEndMS,...
%   plotRTStartMS)
function [T, maxDate] = getTable(dataDirName, tableDataName)

  if ~isfile(tableDataName)
      varNames = {'animal', 'date', 'numStim', ...
          'meanPowerMW', 'maxPowerMW', ...
          'stimCorrects', 'stimFails', 'stimEarlies', ...
          'numNoStim', 'noStimCorrects', 'noStimFails', 'noStimEarlies', ...
          'RTWindowMS', 'pHit', 'pFA', 'dPrime', 'c', ...
          'stimPHit', 'stimPFA', 'stimDPrime', 'stimC', ...
          'noStimPHit', 'noStimPFA', 'noStimDPrime', 'noStimC', ...
          'peakMinMS', 'peakMaxMS', 'rampMS', 'kernelCI', 'kernelPeak', ...
          'hitKernel', 'failKernel', 'earlyKernel', 'RTKernel', 'SRTKernel', 'randomKernel', ...
          'correctRTs', 'failRTs', 'earlyRTs'};
      varTypes = {'string', 'string', 'uint32', ...
          'double', 'double', ...
          'uint32', 'uint32', 'uint32', ...
          'uint32', 'uint32', 'uint32', 'uint32', ...
          'double', 'double', 'double',  'double', 'double', ...
          'double', 'double', 'double',  'double', ...
          'double', 'double', 'double',  'double', ...
          'double', 'double', 'uint32', 'double', 'double', ...
          'cell', 'cell', 'cell', 'cell', 'cell', 'cell', ...
          'cell', 'cell', 'cell'};
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
          if length(fileName) > 10 || str2double(fileName(1:4)) > 2020	% other sorts of data in folder
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
    row.meanPowerMW = 0.0;
    row.maxPowerMW = 0.0;
    row.stimCorrects = 0;
    row.stimFails = 0;
    row.stimEarlies = 0;
    row.numNoStim = 0;
    row.noStimCorrects = 0;
    row.noStimFails = 0;
    row.noStimEarlies = 0;
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
    row.earlyKernel = {0};
    row.randomKernel = {0};
    row.correctRTs = {0};
    row.failRTs = {0};
    row.earlyRTs = {0};
    row.rampMS = 0;                                           % no data for the requested row
    % load the data
    clear file trials dParams
    load([dataDirName row.animal '/MatFiles/' row.date]);     %#ok<LOAD>
    if exist('trials', 'var') && isfield(trials, 'trial')     %#ok<NODEF>
        row.rampMS = trials(1).trial.visualRampDurMS;
        trials = validateTrials(trials);
        row = getKernels(file, trials, row);    
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