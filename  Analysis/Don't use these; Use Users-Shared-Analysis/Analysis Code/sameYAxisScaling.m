function yLimits = sameYAxisScaling(subPlotRows, subPlotColumns, subPlots, prescribedLimits)

% put a set of subplots on the same y-axis scale

  plots = length(subPlots);
  if plots < 1
      return;
  end
  ax = cell(1, plots);
  allYLimits = zeros(plots, 2);
  for p = 1:plots
      subplot(subPlotRows, subPlotColumns, subPlots(p));
      ax{p} = gca;
      allYLimits(p, :) = ylim;
  end
  if nargin >= 4                % if we're given limits, just apply them
    yLimits = prescribedLimits;
  else                          % otherwise, find the max and min across all the plots
    yLimits = [min(allYLimits(:, 1)), max(allYLimits(:, 2))];
  end
  for p = 1:plots
      ylim(ax{p}, yLimits);
  end
end
