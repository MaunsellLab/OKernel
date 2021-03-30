function [respLimitsMS, newIndices, fitCum, endCumTimeMS] = getResponseLimits(file, trials, indices)
%
% The the cumulative RT distribution to find the time range in which responses actually occurred.
% Adjust the contents of indices to exclude RT outside this range.  We fit a logistic function to the FA detrended
% cumulative response function and then take upper and lower points of the response function to set the time when
% the response window starts and ends.
%
% The input indices are logical arrays for the relevant trials within the trials argument.  The output indices
% are also logical arrays, adjusted using the response window determined.
%
  upperLimit = 0.975;   % upper point on cumulative hit function, ending response window
  lowerLimit = 0.025;   % lower point on cumulative hit function, starting response window
  allRTs = [[trials(indices.correct).reactTimeMS], [trials(indices.fail).reactTimeMS],  [trials(indices.early).reactTimeMS]];
  if isempty(allRTs)
    respLimitsMS = [];
    newIndices = [];
    fitCum = [];
    endCumTimeMS = [];
    return
  end
  % Set the time limits, and clip off the early and late RTs, keeping track of their numbers
%   startTime = -1000;
  endCumTimeMS = 2000;
  startTime = -file.preStimMinMS;
  endTime = min(file.responseLimitMS, endCumTimeMS);
  numTotal = length(allRTs);
  earlyRTs = allRTs < startTime;
  lateRTs = allRTs >= endTime;
  RTs = allRTs(~earlyRTs & ~lateRTs);
  % make a cumulative distribution, using the early and late numbers to set the start and end points
  RTDist = zeros(1, endTime - startTime);
  for i = 1:length(RTs)
    bin = RTs(i) - startTime + 1;
    RTDist(bin) = RTDist(bin) + 1;
  end
	numEarly = sum(earlyRTs);
  cumDist = (cumsum(RTDist) + numEarly) / numTotal;
  
  % Fit a line to the early part of the cumulative RT distribution, then deslope the distribution
  b = polyfit(1:file.preStimMinMS, cumDist(1:file.preStimMinMS), 1);
  xData = 1:length(cumDist);
  deSloped = cumDist - xData * b(1) - b(2);
  
  % Fit function: a: logistic infliction point, b: logistic minimum amplitude, logistic maximum amplitude
  % d: logistic Hill's slope. The Hill's slope refers to the steepness (positive or negative) of the curve. 
  opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
  ft =  fittype('c + (b - c)/(1 + (x/a)^d)'); 
  maxY = max(deSloped);
  minY = min(deSloped);
  opts.StartPoint = [(minY + maxY) / 2.0, minY, maxY, 3];
  opts.Lower = [0, -abs((minY + maxY) / 2.0), maxY / 2.0, 0];
  opts.Upper = [length(deSloped), maxY, 2 * maxY, 1000];
  [fitResult, ~] = fit(xData', deSloped', ft, opts ); 

  fitCum = (fitResult.c + (fitResult.b - fitResult.c) ./ (1 + (xData ./ fitResult.a).^fitResult.d)) + ...
      xData * b(1) + b(2);
  respLimitsMS = fitResult.a .* exp(log(1 ./ [upperLimit, lowerLimit] - 1) ./ fitResult.d) + startTime;
  trialEnds = [trials(:).trialEnd];
  RTs = [trials(:).reactTimeMS];
  if length(RTs) > length(trialEnds) 
    RTs = zeros(1, length(trialEnds));
    for t = 1:length(trialEnds)
      RTs(t) = trials(t).reactTimeMS(1);
    end
  end
  newIndices.correct = trialEnds == 0 & RTs >= respLimitsMS(1) & RTs < respLimitsMS(2);
  newIndices.early = trialEnds == 1 & RTs < respLimitsMS(1);
  newIndices.fail = trialEnds == 2 & RTs >= respLimitsMS(2);
end
