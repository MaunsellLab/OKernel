function preProcessAll
% Preprocess OKernel data for analysis.  This should be run whenever key analysis features are changed
% Starting from scratch, prepare all data files for analysis.  We scan through every file doing the following:
%  1) Find the response window to use, and assign trials to hits, misses, earlies
%  2) Add a row for the file in the master table
%  3) Make the stim profiles that are needed for bootstrapping kernels
%
% The data set processed is seleceted using whichData.m.  Currently supported are 'JDC' and 'JJC' for 
% Julian Day-Cooney's and Jackson Cone's data

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
%           [row, stimProfiles] = doOneFile(dataDirName, '1145', '2020-03-21');
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
	if  (strcmp(animalName, '902') && fileName < "2019-09-13") || ...       // exclude training days
      (strcmp(animalName, '905') && fileName < "2019-09-24") || ...
      (strcmp(animalName, '1112') && fileName < "2020-02-04") || ...
      (strcmp(animalName, '1145') && fileName < "2020-02-10") || ...
      (strcmp(animalName, '1150') && fileName < "2020-01-16") || ...
      (strcmp(animalName, '1218') && fileName < "2020-03-23") || ...
      (strcmp(animalName, '1220') && fileName < "2020-04-07") || ...
      (strcmp(animalName, '1223') && fileName < "2020-02-16") || ...
      (strcmp(animalName, '1257') && fileName < "2020-04-04")      
    row = [];
    stimProfiles = [];
    return;
  end
  load([dataDirName animalName '/MatFiles/' fileName]); %#ok<LOAD>
  if ~exist('trials', 'var') || ~isfield(trials, 'trial') || ~(str2double(animalName) == file.subjectNumber) %#ok<NODEF>
    row = [];
    stimProfiles = [];
    return;
  end
	trials = validateTrials(trials);                        % check trialEnds, etc.
  row = initializeRow(animalName, fileName, trials(1).trial.visualRampDurMS);
	[row, stimProfiles] = getKernels(file, trials, row);      % compile the kernels for this file
end

%%
function animalNames = getAnimalNames()

    [dirName, ~, projName] = whichData();
    dirStructs = dir(dirName);                              % data directory contents
    animalNames = {dirStructs(:).name};                       % data directory file names
    numFiles = length(animalNames);
    validFiles = false(1, numFiles);                        % find only files with number names (animals)
    for f = 1:numFiles
      validFiles(f) = isstrprop(animalNames{f}(1), 'digit') && length(animalNames{f}) > 1;
    end
    animalNames = animalNames(validFiles);
    fileValues = str2double(animalNames);                   % sort the files numerically
    switch projName                                         % looking for only certain animals
      case 'JDC'
        validFiles = fileValues < 1400;
      case 'JJC'
        validFiles = fileValues >= 1400;
      otherwise
        fprintf('preProcessAll:getAnimalNames: unrecognized project name: %s', projName);
    end
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
  [stimIndices, trials] = getStimIndices(trials);
  if sum(stimIndices) == 0
    row = [];
    stimProfiles = [];
    return
  end
  trialStructs = [trials(:).trial];
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
    row = []; stimProfiles = [];
    return;
  end
  % find the response interval and get a modified set of indices that limits hits to only that interval
	[respLimitsMS, theIndices, ~, ~] = getResponseLimits(file, trials, theIndices);
  row.startRT = respLimitsMS(1);
  row.endRT = respLimitsMS(2);
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
  if abs(row.dPrime) == Inf
    row.dPrime = NaN; 
  end

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
    if abs(row.noStimDPrime) == Inf
      row.noStimDPrime = NaN; 
    end
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

  % Get the early kernel.  Eliminate trials that start before the end of the ramping stimulus.
	earlyIndices = indices.early & (preStimMS + RTs + plotRTStartMS > 200);
  if row.stimEarlies > 0
      profiles = getStimProfiles(trials(earlyIndices), plotRTStartMS, plotRTStartMS + plotEndMS - plotStartMS, true, true);
      earlySum(:) = sum(profiles, 1);
      row.stimEarlies = size(profiles, 1);                  % getStimProfiles might reject some trials as too short
      row.earlyKernel = {earlySum / row.stimEarlies};
  end
  
  % If we have corrects and fails, save the hit, fail, total and random kernels.
  row.kernelCI = sqrt(hitCI^2 + failCI^2);
  if ~isempty(hitKernel) && ~isempty(failKernel)
    row.kernelPeak =  max(abs((hitKernel(1:plotEndMS - plotStartMS) - failKernel(1:plotEndMS - plotStartMS))) / row.kernelCI);
    row.randomKernel = getRandomKernel(row.stimCorrects, row.stimFails, trialStructs(1).pulseDurMS, plotEndMS - plotStartMS);
  end

  % add to the RT distributions
  row.correctRTs = {[trials(theIndices.correct).reactTimeMS]};
  row.earlyRTs = {[trials(theIndices.early).reactTimeMS]};
  failRTs = [trials(theIndices.fail).reactTimeMS];
  failRTs(failRTs < 0 | failRTs > 100000) = 100000;     % include fails in count, but don't let them display on plot
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
  row.correctRTs = {0};
  row.failRTs = {0};
  row.earlyRTs = {0};
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
    'numStim', 'stimCorrects', 'stimFails', 'stimEarlies', ...
    'numNoStim', 'noStimCorrects', 'noStimFails', 'noStimEarlies', ...
    'pHit', 'pFA', 'dPrime', 'c', ...
    'stimPHit', 'stimPFA', 'stimDPrime', 'stimC', 'noStimPHit', 'noStimPFA', 'noStimDPrime', 'noStimC', ...
    'kernelCI', 'kernelPeak', 'hitKernel', 'failKernel', 'earlyKernel', 'RTKernel', 'SRTKernel', 'randomKernel', ...
    'startRT', 'endRT', 'correctRTs', 'failRTs', 'earlyRTs'};
  types = {'string', 'string', 'uint32', 'double', 'double', 'double', ...
    'uint32', 'uint32', 'uint32', 'uint32', ...
    'uint32', 'uint32', 'uint32', 'uint32', ...
    'double', 'double',  'double', 'double', ...
    'double', 'double', 'double',  'double', 'double', 'double', 'double',  'double', ...
    'double', 'double', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', ...
    'uint32', 'uint32', 'cell', 'cell', 'cell'};
end

 