function xLimits = sameXAxisScaling(subPlotRows, subPlotColumns, subPlots, prescribedLimits)

% put a set of subplots on the same y-axis scale

  plots = length(subPlots);
  if plots < 1
      return;
  end
  ax = cell(1, plots);
  allXLimits = zeros(plots, 2);
  for p = 1:plots
      subplot(subPlotRows, subPlotColumns, subPlots(p));
      ax{p} = gca;
      allXLimits(p, :) = xlim;
  end
  if nargin >= 4                % if we're given limits, just apply them
    if prescribedLimits(1) == inf
      prescribedLimits(1) = min(allXLimits(:, 1));
    end
    if prescribedLimits(2) == inf
      prescribedLimits(2) = max(allXLimits(:, 2));
    end
    xLimits = prescribedLimits;
  else                          % otherwise, find the max and min across all the plots
    xLimits = [min(allXLimits(:, 1)), max(allXLimits(:, 2))];
  end
  for p = 1:plots
      xlim(ax{p}, xLimits);
  end
end
