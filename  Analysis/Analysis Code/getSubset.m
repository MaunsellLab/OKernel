function [U, maxDate] = getSubset(mode, dataDirName, tableDataName, limits)
  % Return a subset of files based on selection criteria.
  % If necessary, a new table of summary results will be compiled.
    
  [T, maxDate] = getTable(dataDirName, tableDataName);                % get table of all sessions
  switch mode
    case 'normal'
      controls = controlSessions(T);                                  % logic array flagging control sesions
      valid = any(T.rampMS == limits.rampMS & T.kernelCI > 0, 2) & ~controls; % empty entries have zero for kernelCI
    case 'control'
       controls = controlTestSessions(T);                         	  % logic array flagging control test sesions
       valid = any(T.rampMS == limits.rampMS & T.kernelCI > 0, 2) & controls; % empty entries have zero for kernelCI
   case 'prePostControl'
       controls = prePostControlSessions(T);                         	% logic array flagging control test sesions
       valid = any(T.rampMS == limits.rampMS & T.kernelCI > 0, 2) & controls; % empty entries have zero for kernelCI
    otherwise
      fprintf('getSubset: unrecognized mode');
  end
  if ~strcmp(limits.animal, 'All')
    valid = valid & T.animal == limits.animal;
  end
  if ~isempty(limits.oneDay)
    valid = valid & T.date == limits.oneDay;
  end
  valid = valid & (T.stimCorrects > limits.minTrials);                % enough correct trials
  valid = valid & T.stimFails > limits.minTrials;                    	% enough fail trials
  valid = valid & T.kernelPeak > limits.criterion;
  
  % d' decrement limits
  if limits.minDec ~= -1
    valid = valid & double(T.noStimDPrime) - double(T.stimDPrime) >= limits.minDec;     % min behavioral decrement in hit rate
  end
  U = T(valid, :);
  U = withoutTrainingDays(U);                                         % remove sessions from training periods
  
  % remove any animals with too few sessions for any ramp duration
  if isempty(limits.oneDay) && limits.minSessions > 0
    animals = unique(U.animal);
    rampDurs = unique(limits.rampMS);
    for r = 1:length(rampDurs)
      for a = 1:length(animals)
        animalRows = U.animal == animals{a} & U.rampMS == rampDurs(r);
        if sum(animalRows) < limits.minSessions
          U = U(~animalRows, :);
        end
      end
    end
  end
  
  if size(U, 1) == 0
    if length(limits.rampMS) == 1
      fprintf('getSubset: No valid sessions found for ''%s'' for animal ''%s'' on %d ms ramps, minTrials %d and decrement %.2f\n', ...
        mode, limits.animal{1}, limits.rampMS, limits.minTrials, limits.minDec);
    else
      fprintf('getSubset: No valid sessions found for ''%s'' for animal ''%s'' on multiple ramps, minTrials %d and decrement %.2f\n', ...
        mode, limits.animal{1}, limits.minTrials, limits.minDec);
    end
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

% function indices = getAllIndices(trials, eotCodes, stimIndices)
% %
% % Return logical arrays specifying valid indices for hit, fail and early trials.  The lists are modified
% % according to trialType, which can select all, stimulated, or unstimulated trials.
% 
%   trialStructs = [trials(:).trial];                           % extract trialStructs to an array for access
%   RTs = [trials(:).reactTimeMS];                              % get all trial RTs
%   selectIndices = true(1, length(stimIndices));            % take all trials
% 	indices.correct = selectIndices & eotCodes == 0;            % get hit trials
%   indices.early = selectIndices & eotCodes == 1;              % get early trials
%   indices.fail = selectIndices & eotCodes == 2;               % get fail trials
%   
%   % for earlies, we must eliminate those that don't have enough time before the FA to avoid sampling before 
%   % the end of the ramp (which is fixed at 200 ms). 
%   preStimMS = [trialStructs(:).preStimMS];                  	% get preStim times for each trial
%   [~, ~, plotRTStartMS] = plotLimits();                       % get the limits for the plots we will make
%   indices.early = (selectIndices & eotCodes == 1) & (preStimMS + RTs + plotRTStartMS > 200);
% end

%%
function row = getKernels(file, trials, row)
%
% Compute kernels for one session
%
  trialStructs = [trials(:).trial];
  meanPower = [trials(:).meanPowerMW];                        % get power applied for each trial                        
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
      eotCodes(t) = trials(t).trialEnd;
  end
  RTs = [trials(:).reactTimeMS];                              % get all trial RTs  
	preStimMS = [trialStructs(:).preStimMS];                  	% get preStim times for each trial  

  % calculate the overall d'
  allIndices = allIndices(trials, eotCodes, stimIndices);
  if sum(allIndices.correct | allIndices.fail | allIndices.early) < 10
    return;
  end
  % find the response interval and get a modified set of indices that limits hits to only that interval
	[respLimitsMS, allIndices, ~, ~] = getResponseLimits(file, trials, allIndices);
  row.RTWindowMS = diff(respLimitsMS);
  row.noStimCorrects = sum(allIndices.correct & ~stimIndices);
  row.noStimFails = sum(allIndices.fail & ~stimIndices);
  row.noStimEarlies = sum(allIndices.early & ~stimIndices); 
  row.numNoStim = row.noStimCorrects + row.noStimFails + row.noStimEarlies;  % don't count kEOTIgnored
  row.stimCorrects = sum(allIndices.correct & stimIndices);
  row.stimFails = sum(allIndices.fail & stimIndices);
  row.stimEarlies = sum(allIndices.early & stimIndices);
  row.numStim = row.stimCorrects + row.stimFails + row.stimEarlies;                  % don't count kEOTIgnored

  % find the performance across stim and nostim trials combined
  hitRate = sum(allIndices.correct) / (sum(allIndices.correct) + sum(allIndices.fail));
  rateEarly = earlyRate(file, trials, allIndices.correct, allIndices.fail, allIndices.early);
  row.pFA = 1.0 - exp(-rateEarly * row.RTWindowMS / 1000.0);
  row.pHit = (hitRate - row.pFA) / (1.0 - row.pFA);
  [row.dPrime, row.c] = dprime(row.pHit, row.pFA);
  
  % calculate the nostim trial d' using the all-trial pFA and the all-trial response window
  indices.correct = allIndices.correct & ~stimIndices;
  indices.fail = allIndices.fail & ~stimIndices;
  indices.early = allIndices.early & ~stimIndices;
  if sum(indices.correct | indices.fail | indices.early) > 0
    hitRate = sum(indices.correct) / (sum(indices.correct) + sum(indices.fail));
    rateEarly = earlyRate(file, trials, indices.correct, indices.fail, indices.early);
    row.noStimPFA = 1.0 - exp(-rateEarly * row.RTWindowMS / 1000.0);
    row.noStimPHit = (hitRate - row.pFA) / (1.0 - row.pFA);         % using overall pFA
    [row.noStimDPrime, row.noStimC] = dprime(row.noStimPHit, row.pFA);
%     fprintf('  No Stim: pFA = %.2f, pHit = %.2f d-prime = %.1f, c = %.1f\n', ...
%               row.noStimPFA, row.noStimPHit, row.noStimDPrime, row.noStimC);
  end
    
  % calculate the stim trial d' using the all-trial pFA and the all-trial response window
  indices.correct = allIndices.correct & stimIndices;
  indices.fail = allIndices.fail & stimIndices;
  indices.early = allIndices.early & stimIndices;
  if sum(indices.correct | indices.fail | indices.early) > 0
    hitRate = sum(indices.correct) / (sum(indices.correct) + sum(indices.fail));
    rateEarly = earlyRate(file, trials, indices.correct, indices.fail, indices.early);
    row.stimPFA = 1.0 - exp(-rateEarly * row.RTWindowMS / 1000.0);
    row.stimPHit = (hitRate - row.pFA) / (1.0 - row.pFA);   % using overall pFA
    [row.stimDPrime, row.stimC] = dprime(row.stimPHit, row.pFA);
%     fprintf('  Stim:    pFA = %.2f, pHit = %.2f d-prime = %.1f, c = %.1f  delta: %.5f\n', ...
%               row.stimPFA, row.stimPHit, row.stimDPrime, row.stimC, row.noStimDPrime - row.stimDPrime);
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
  row.correctRTs = {[trials(allIndices.correct).reactTimeMS]};
  row.earlyRTs = {[trials(allIndices.early).reactTimeMS]};
  failRTs = [trials(allIndices.fail).reactTimeMS];
  failRTs(failRTs < 0 | failRTs > 100000) = 100000;     % include fails in count, but don't let them display on plot
  row.failRTs = {failRTs};
end

%% getTable -- load or create the table with all the summary data
% function [T, maxDate] = getTable(dataDirName, tableDataName, peakMinMS, peakMaxMS, plotStartMS, plotEndMS,...
%   plotRTStartMS)
function [T, maxDate] = getTable(dataDirName, tableDataName)

  if ~isfile(tableDataName)
      varNames = {'animal', 'date', 'numStim', 'stimCorrects', 'stimFails', 'stimEarlies', ...
          'numNoStim', 'noStimCorrects', 'noStimFails', 'noStimEarlies', ...
          'RTWindowMS', 'pHit', 'pFA', 'dPrime', 'c', ...
          'stimPHit', 'stimPFA', 'stimDPrime', 'stimC', ...
          'noStimPHit', 'noStimPFA', 'noStimDPrime', 'noStimC', ...
          'peakMinMS', 'peakMaxMS', 'rampMS', 'kernelCI', 'kernelPeak', ...
          'hitKernel', 'failKernel', 'earlyKernel', 'RTKernel', 'SRTKernel', 'randomKernel', ...
          'correctRTs', 'failRTs', 'earlyRTs'};
      varTypes = {'string', 'string', 'uint32', 'uint32', 'uint32', 'uint32', ...
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
        trials = validateTrials(row, trials);
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
