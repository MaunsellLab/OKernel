function preProcessAll
% Preprocess SCernel data for analysis.  This should be run whenever key analysis features are changed
% Starting from scratch, prepare all data files for analysis.  We scan through every file doing the following:
%  1) Find the response window to use, and assign trials to hits, misses, earlies
%  2) Add a row for the file in the master table
%  3) Make the stim profiles that are needed for bootstrapping kernels
%
% The data set processed is selected using whichData.m. 

  [dataDirName, tableName] = whichData();
  [varNames, varTypes] = tableNamesAndTypes();
  T = table('size', [0, length(varNames)], 'variableTypes', varTypes, 'variableNames', varNames);
  % count the number of sessions for online display
  animalNames = getAnimalNames();
  numAnimals = length(animalNames);
  numSessions = 0;
  for a = 1:numAnimals
      fileNames = getFileNames(animalNames{a});
      numSessions = numSessions + length(fileNames);
  end
  % compile each session
  session = 1;
  for a = 1:numAnimals
      fileNames = getFileNames(animalNames{a});
      numFiles = length(fileNames);
      for f = 1:numFiles                                  	% for each file
          [~, fileName, ~] = fileparts(fileNames{f});      	% get file name
          fprintf('Session %4d of %4d: %4s %s\n', session, numSessions, animalNames{a}, fileName);
          session = session + 1;
          if length(fileName) > 10
             continue; 
          end
          [row, stimProfiles] = doOneFile(dataDirName, animalNames{a}, fileName);
          if ~isempty(row)
            T = [T; row];   %#ok<AGROW>
          end
          if ~isempty(stimProfiles)
            if ~(isfolder([dataDirName, ' Analysis/Mat Files/Stim Profiles/']))
              mkdir([dataDirName, ' Analysis/Mat Files/Stim Profiles/'])
            end
            if ~(isfolder([dataDirName, ' Analysis/Mat Files/Stim Profiles/', animalNames{a}]))
              mkdir([dataDirName, ' Analysis/Mat Files/Stim Profiles/', animalNames{a}])
            end
            save([dataDirName, ' Analysis/Mat Files/Stim Profiles/', animalNames{a}, '/', fileName], 'stimProfiles'); 
          end
      end
  end
  save([dataDirName, tableName], 'T'); 
end

%%
function [row, stimProfiles] = doOneFile(dataDirName, animalName, fileName)
%
  load([dataDirName animalName '/MatFiles/' fileName]); %#ok<LOAD>
  if ~exist('trials', 'var') || ~isfield(trials, 'trial') || ~(str2double(animalName) == file.subjectNumber) %#ok<NODEF>
    row = [];
    stimProfiles = [];
    return;
  end
	trials = validateTrials(trials);                          % check trialEnds, etc.
  row = initializeRow(animalName, fileName, trials(1).trial.visualRampDurMS);
	[row, stimProfiles] = getKernels(file, trials, row);      % compile the kernels for this file
end

%%
function animalNames = getAnimalNames()

    [dirName] = whichData();
    dirStructs = dir(dirName);                              % data directory contents
    animalNames = {dirStructs(:).name};                       % data directory file names
    numFiles = length(animalNames);
    validFiles = false(1, numFiles);                        % find only files with number names (animals)
    for f = 1:numFiles
      validFiles(f) = isstrprop(animalNames{f}(1), 'digit') && length(animalNames{f}) > 1;
    end
    animalNames = animalNames(validFiles);
    fileValues = str2double(animalNames);                   % sort the files numerically
    validFiles = fileValues >= 1400;
    animalNames = animalNames(validFiles);
    fileValues = fileValues(validFiles);
    [~, indices] = sort(fileValues);
    animalNames = animalNames(indices);
end

%%
function fileNames = getFileNames(animalName)

    dataDirName = whichData();
    dirStructs = dir([dataDirName animalName '/MatFiles/']);	% data directory contents
    fileNames = {dirStructs(:).name};                             % data directory file names
    numFiles = length(fileNames);
    validFiles = false(1, numFiles);                              % find only files with number names (animals)
    for f = 1:numFiles
      validFiles(f) = isstrprop(fileNames{f}(1), 'digit') && length(fileNames{f}) > 1;
      validFiles(f) = validFiles(f) & length(fileNames{f}) == 14;
    end
    fileNames = fileNames(validFiles);
    fileNames = sort(fileNames);
end

%%
function [row, stimProfiles] = getKernels(file, trials, row)
%
% Compute kernels for one session
%
  trialStructs = [trials(:).trial];
  meanPower = [trials(:).meanPowerMW];                        % get power applied for each trial 
  row.meanPowerMW = mean(meanPower);
  row.maxPowerMW = max(meanPower);
  stimIndices = meanPower > 0;                                % trials with opto stimulation
  if isfield(trialStructs, 'pulseContrast')                   % get rid of any trials with reduced opto power
       stimIndices = stimIndices & [trialStructs.pulseContrast] == 1;	% but only those with contrast == 1
  end
  if sum(stimIndices) == 0                                    % no stimulation with full opto contrast
    row = []; stimProfiles = [];
    return;
  end
  
  % we only consider the range between the first and last stimulated trials
  firstStimIndex = find(stimIndices > 0, 1);                 	% first stimulated trial
	lastStimIndex = find(stimIndices > 0, 1, 'last');          	% last stimulated trial
  trials = trials(firstStimIndex:lastStimIndex);
  trialStructs = trialStructs(firstStimIndex:lastStimIndex);
  stimIndices = stimIndices(firstStimIndex:lastStimIndex);
  
  eotCodes = [trials(:).trialEnd];                            % extract all trialEnds
  RTs = [trials(:).reactTimeMS];                              % extract all trial RTs  
	preStimMS = [trialStructs(:).preStimMS];                    % extract preStim times for each trial  

  % calculate the overall d'
  theIndices = allIndices(trials, eotCodes);
  if sum(theIndices.correct | theIndices.fail | theIndices.early) < 10
    row = []; stimProfiles = [];
    return;
  end
    
  % Get Visual Stimulus Levels to subset data by top-up versus low (tested) contrast
  contrasts = [trialStructs(:).visualStimValue];
  testContIdx = contrasts == min(contrasts); % We currently test the min contrast presented
  topUpContIdx = contrasts == max(contrasts); % We present high contrasts just to keep em happy
  
  % find the response interval and get a modified set of indices that limits hits to only that interval
  [respLimitsMS, theIndices, ~, ~] = getResponseLimits(file, trials, theIndices);
  row.startRT = respLimitsMS(1);
  row.endRT = respLimitsMS(2);
  row.RTWindowMS = diff(respLimitsMS);
  
  % Basic metrics for TopUp Contrast Trials
  row.topUpCorrects = sum(theIndices.correct & ~stimIndices & topUpContIdx);
  row.topUpFails = sum(theIndices.fail & ~stimIndices & topUpContIdx);
  row.topUpEarlies = sum(theIndices.early & ~stimIndices & topUpContIdx);  % This is for all trials but thats fine bc no vis stim occurs
  row.numTopUp = row.topUpCorrects + row.topUpFails + row.topUpEarlies;
  
  % Basic metrics for no stim trials @ tested contrast
  row.noStimCorrects = sum(theIndices.correct & ~stimIndices & testContIdx);
  row.noStimFails = sum(theIndices.fail & ~stimIndices & testContIdx);
  row.noStimEarlies = sum(theIndices.early & ~stimIndices & testContIdx);  % This is for all trials but thats fine bc no vis stim occurs
  row.numNoStim = row.noStimCorrects + row.noStimFails + row.noStimEarlies;  % don't count kEOTIgnored
  
  % Basic metric for stim trials @ tested contrast
  row.stimCorrects = sum(theIndices.correct & stimIndices & testContIdx);
  row.stimFails = sum(theIndices.fail & stimIndices & testContIdx);
  row.stimEarlies = sum(theIndices.early & stimIndices & testContIdx);
  row.numStim = row.stimCorrects + row.stimFails + row.stimEarlies;                  % don't count kEOTIgnored

  % Overall Earlies and Early Rate Across All trials
  rateEarly = earlyRate(file, trials, theIndices.correct, theIndices.fail, theIndices.early);
  row.pFA = 1.0 - exp(-rateEarly * row.RTWindowMS / 1000.0);
  
  % find the performance on Top-Up trials 
  indices.correct = theIndices.correct & ~stimIndices & topUpContIdx;
  indices.fail = theIndices.fail & ~stimIndices & topUpContIdx;
  topUpHitRate = sum(indices.correct) / (sum(indices.correct) + sum(indices.fail));
  row.topUpPHit = (topUpHitRate - row.pFA) / (1.0 - row.pFA);
  
  [row.topUpDPrime, row.topUpC] = dprime(row.topUpPHit, row.pFA, true);
  if abs(row.topUpDPrime) == Inf
    row.row.topUpDPrime = NaN; 
  end

  % calculate the nostim trial d' using the all trial pFA & the all-trial response window
  indices.correct = theIndices.correct & ~stimIndices & testContIdx;
  indices.fail = theIndices.fail & ~stimIndices & testContIdx;
  indices.early = theIndices.early & ~stimIndices & testContIdx;
  
  if sum(indices.correct | indices.fail | indices.early) > 0
    hitRate = sum(indices.correct) / (sum(indices.correct) + sum(indices.fail));
    rateEarly = earlyRate(file, trials, indices.correct, indices.fail, indices.early);
    row.noStimPFA = 1.0 - exp(-rateEarly * row.RTWindowMS / 1000.0);
    row.noStimPHit = (hitRate - row.pFA) / (1.0 - row.pFA);         % using overall pFA
    [row.noStimDPrime, row.noStimC] = dprime(row.noStimPHit, row.pFA, true);
    if abs(row.noStimDPrime) == Inf
      row.noStimDPrime = NaN; 
    end
  end
    
  % calculate the stim trial d' using the all-trial pFA and the all-trial response window
  indices.correct = theIndices.correct & stimIndices & testContIdx;
  indices.fail = theIndices.fail & stimIndices & testContIdx;
  indices.early = theIndices.early & stimIndices & testContIdx;
  
  if sum(indices.correct | indices.fail | indices.early) > 0
    hitRate = sum(indices.correct) / (sum(indices.correct) + sum(indices.fail));
    rateEarly = earlyRate(file, trials, indices.correct, indices.fail, indices.early);
    row.stimPFA = 1.0 - exp(-rateEarly * row.RTWindowMS / 1000.0);
    row.stimPHit = (hitRate - row.pFA) / (1.0 - row.pFA);   % using overall pFA
    [row.stimDPrime, row.stimC] = dprime(row.stimPHit, row.pFA, true);
    if abs(row.stimDPrime) == Inf
      row.stimDPrime = NaN; 
    end
  end
  
  % get the various kernels.  We use the indices set up in the previous block that include only stimulated trials.
  [plotStartMS, plotEndMS, plotRTStartMS] = plotLimits();     % get the limits for the plots we will display  
  
  % get the hit kernel.  We use the indices set up to have only stimulated trials within the detected response window.
  profiles = getStimProfiles(trials(indices.correct), plotStartMS, plotEndMS, true, false);
  hitSums(:) = sum(profiles, 1);
  hitKernel = hitSums / row.stimCorrects;
  row.hitKernel = {hitKernel};
  hitCI = stimCI(row.stimCorrects);

  % get the RT aligned kernel
  profiles = getStimProfiles(trials(indices.correct), plotRTStartMS, plotRTStartMS + plotEndMS - plotStartMS, true, true);
  RTSum(:) = sum(profiles, 1);
  row.RTKernel = {RTSum / row.stimCorrects};

  % get the Stim & RT aligned hit kernel
  profiles = getStimRTProfiles(trials(indices.correct), plotStartMS, plotEndMS);
  SRTSum(:) = sum(profiles, 1);
  row.SRTKernel = {SRTSum / row.stimCorrects};

  % get the fail kernel
  profiles = getStimProfiles(trials(indices.fail), plotStartMS, plotEndMS, true, false);
  failSum(:) = sum(profiles, 1);
  failKernel = failSum / row.stimFails;
  row.failKernel = {failKernel};
  failCI = stimCI(row.stimFails);

  % Get the early kernel.  
  if row.stimEarlies > 0
      profiles = getStimProfiles(trials(indices.early), plotRTStartMS, plotRTStartMS + plotEndMS - plotStartMS, true, true);
      earlySum(:) = sum(profiles, 1);
      row.stimEarlies = size(profiles, 1);                  % getStimProfiles might reject some trials as too short
      row.earlyKernel = {earlySum / row.stimEarlies};
  end
  
  % If we have corrects and fails, save the hit, fail, total and random kernels.
  row.kernelCI = sqrt(hitCI^2 + failCI^2);
  if ~isempty(hitKernel) && ~isempty(failKernel)
    row.kernelPeak =  max(abs((hitKernel(1:plotEndMS - plotStartMS) - failKernel(1:plotEndMS - plotStartMS))) / row.kernelCI);
    row.randomKernel = {getRandomKernel(row.stimCorrects, row.stimFails, trialStructs(1).pulseDurMS, plotEndMS - plotStartMS)};
  end

  % add to the RT distributions (Do this for each stim type, can be re-combined later if needed)
  row.topUpCorrectRTs = {[trials(theIndices.correct & ~stimIndices & topUpContIdx).reactTimeMS]};
  row.stimCorrectRTs = {[trials(theIndices.correct & stimIndices & testContIdx).reactTimeMS]};
  row.noStimCorrectRTs = {[trials(theIndices.correct & ~stimIndices & testContIdx).reactTimeMS]};
  
  row.topUpEarlyRTs = {[trials(theIndices.early & ~stimIndices & topUpContIdx).reactTimeMS]};
  row.stimEarlyRTs = {[trials(theIndices.early & stimIndices & testContIdx).reactTimeMS]};
  row.noStimEarlyRTs = {[trials(theIndices.early & ~stimIndices & testContIdx).reactTimeMS]};
  
  failRTs = [trials(theIndices.fail).reactTimeMS];
  failRTs(failRTs < 0 | failRTs >= 10000) = 100000;     % include fails in count, but don't let them display on plot
  row.failRTs = {failRTs};
    
  % get the hit profiles
  hitIndices = theIndices.correct & stimIndices;
  hitProfiles = getStimProfiles(trials(hitIndices), plotStartMS, plotEndMS, true, false);
  stimProfiles.hitProfiles = normProfiles(hitProfiles);

  % get the RT aligned profiles
  RTProfiles = getStimProfiles(trials(hitIndices), plotRTStartMS, plotRTStartMS + plotEndMS - plotStartMS, true, true);
  stimProfiles.RTProfiles = normProfiles(RTProfiles);

  % get the Stim & RT aligned hit profiles
  stimRTProfiles = getStimRTProfiles(trials(hitIndices), plotStartMS, plotEndMS);
  stimProfiles.stimRTProfiles = normProfiles(stimRTProfiles);

  % get the miss profiles
  missProfiles = getStimProfiles(trials(theIndices.fail & stimIndices), plotStartMS, plotEndMS, true, false);
  stimProfiles.missProfiles = normProfiles(missProfiles);

  % Get the early profiles
  earlyProfiles = getStimProfiles(trials(theIndices.early & stimIndices), plotRTStartMS, plotRTStartMS + plotEndMS - plotStartMS, true, true);
  stimProfiles.earlyProfiles = normProfiles(earlyProfiles);  
end

%%
function row = initializeRow(animalName, fileName, rampMS)
%
  [varNames, varTypes] = tableNamesAndTypes();
  row = table('size', [1, length(varNames)], 'variableTypes', varTypes, 'variableNames', varNames);
  row.animal = animalName;
  row.date = fileName;
	row.rampMS = rampMS;
  row.meanPowerMW = 0.0;
  row.maxPowerMW = 0.0;
  row.RTWindowMS = 0.0;
  row.numTopUp = 0;
  row.topUpCorrects = 0;
  row.topUpFails = 0;
  row.topUpEarlies = 0;
  row.numStim = 0;
  row.stimCorrects = 0;
  row.stimFails = 0;
  row.stimEarlies = 0;
  row.numNoStim = 0;
  row.noStimCorrects = 0;
  row.noStimFails = 0;
  row.noStimEarlies = 0;
  row.pFA = 0.0;
  row.topUpPHit = 0.0;
  row.topUpDPrime = 0.0;
  row.topUpC = 0.0;
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
  row.startRT = 0;
  row.endRT = 0;
  row.topUpCorrectRTs = {0};
  row.stimCorrectRTs = {0};
  row.noStimCorrectRTs = {0};
  row.topUpEarlyRTs = {0};
  row.stimEarlyRTs = {0};
  row.noStimEarlyRTs = {0};
  row.failRTs = {0};
end

%%
function profiles = normProfiles(profiles)

meanPower = (max(max(profiles)) + min(min(profiles))) / 2.0;
profiles(profiles < meanPower) = -1.0;
profiles(profiles >= meanPower) = 1.0;
end

%%
function [names, types] = tableNamesAndTypes()
%
  names = {'animal', 'date', 'rampMS', 'meanPowerMW', 'maxPowerMW', 'RTWindowMS', ...
    'numTopUp', 'topUpCorrects', 'topUpFails', 'topUpEarlies', ...
    'numStim', 'stimCorrects', 'stimFails', 'stimEarlies', ...
    'numNoStim', 'noStimCorrects', 'noStimFails', 'noStimEarlies', ...
    'pFA', 'topUpPHit', 'topUpDPrime', 'topUpC', ...
    'stimPHit', 'stimPFA', 'stimDPrime', 'stimC', 'noStimPHit', 'noStimPFA', 'noStimDPrime', 'noStimC', ...
    'kernelCI', 'kernelPeak', 'hitKernel', 'failKernel', 'earlyKernel', 'RTKernel', 'SRTKernel', 'randomKernel', ...
    'startRT', 'endRT', ...
    'topUpCorrectRTs', 'stimCorrectRTs','noStimCorrectRTs' 'topUpEarlyRTs', 'stimEarlyRTs', 'noStimEarlyRTs', 'failRTs'};
  types = {'string', 'string', 'uint32', 'double', 'double', 'double', ...
    'uint32', 'uint32', 'uint32', 'uint32', ...
    'uint32', 'uint32', 'uint32', 'uint32', ...
    'uint32', 'uint32', 'uint32', 'uint32', ...
    'double', 'double',  'double', 'double', ...
    'double', 'double', 'double',  'double', 'double', 'double', 'double',  'double', ...
    'double', 'double', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', ...
    'uint32', 'uint32', ...
    'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell'};
end

 