function limits = setLimits(theSubset)
  % Set up the default selection criteria
  limits.minSessions = 3;                	% require at least n sessions for each animal
  limits.minTrials = 0;
  limits.criterion = 0;
  limits.minDec = -1;
  limits.minDPrime = -1;
  limits.minAvgDeltaDPrime = 0.1; % Use this to control stringency for data selection
  limits.minSessionDeltaDPrime = 0.1; % Single Session Delta d-prime must exceed this for inclusion 
  limits.maxMeanPowerMW = 0.5;
  limits.animal = {'All'};
  limits.oneDay = [];
  switch theSubset
    case {'Unfiltered', 'unfiltered'}
      limits.minSessions = 0;
      limits.minDec = -1;
      limits.minDPrime = -1;
      limits.minAvgDeltaDPrime = -1;
      limits.maxMeanPowerMW = 1000;
      limits.rampMS = [0, 500];
    case {'All', 'all', 'all steps', 'All Steps', 'All steps', 'all ramps', 'All ramps', 'All Ramps'}
      switch theSubset
        case {'All', 'all'}
          limits.rampMS = [0, 500];
        case {'all steps', 'All Steps', 'All steps'}
          limits.rampMS = 0;
        case {'all ramps', 'All ramps', 'All Ramps'}
          limits.rampMS = 500;
      end
    case {'Example', 'example'}
      limits.rampMS = 0;
      limits.animal = {'902'};
      limits.oneDay = '2019-10-10';
    case {'OneOff', 'oneoff', 'oneOff'}
      limits.rampMS = 0;
      limits.animal = {'1462', '1463'};
    otherwise
      fprintf('getSessionTable: unrecognized table type ''%s''\n', theSubset);
      limits = [];
      return;
  end
end
