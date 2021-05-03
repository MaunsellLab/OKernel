function limits = sameAxisScaling(theAxis, subPlotRows, subPlotColumns, subPlots, prescribedLimits)
%
% put a set of subplots on the same scale
%
  if nargin < 5
    prescribedLimits = [];
  end
  switch theAxis
    case 'x'
      limits = doScaling(@xlim, subPlotRows, subPlotColumns, subPlots, prescribedLimits);
    case 'y'
      limits = doScaling(@ylim, subPlotRows, subPlotColumns, subPlots, prescribedLimits);
    case {'xy', 'yx', 'both'}
      limits = doScaling(@xlim, subPlotRows, subPlotColumns, subPlots, prescribedLimits);
      limits = doScaling(@ylim, subPlotRows, subPlotColumns, subPlots, prescribedLimits);
  end
end

function newLimits = doScaling(limitFunc, subPlotRows, subPlotColumns, subPlots, prescribedLimits)
%
  plots = length(subPlots);
  if plots < 1
      return;
  end
  ax = cell(1, plots);
  allLimits = zeros(plots, 2);
  for p = 1:plots
      subplot(subPlotRows, subPlotColumns, subPlots(p));
      ax{p} = gca;
      allLimits(p, :) = limitFunc();
  end
  if ~isempty(prescribedLimits)                % if we're given limits, just apply them
    if prescribedLimits(1) == inf
      prescribedLimits(1) = min(allLimits(:, 1));
    end
    if prescribedLimits(2) == inf
      prescribedLimits(2) = max(allLimits(:, 2));
    end
    newLimits = prescribedLimits;
  else                          % otherwise, find the max and min across all the plots
    newLimits = [min(allLimits(:, 1)), max(allLimits(:, 2))];
  end
  for p = 1:plots
      limitFunc(ax{p}, newLimits);
  end
end
