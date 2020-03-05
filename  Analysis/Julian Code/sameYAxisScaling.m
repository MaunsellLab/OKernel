function sameYAxisScaling(subPlotRows, subPlotColumns, subPlots) 
% put a set of subplots on the same y-axis scale

    plots = length(subPlots);
    if plots < 2
        return;
    end
    ax = cell(1, plots);
    yLimits = zeros(plots, 2);
    for p = 1:plots
        subplot(subPlotRows, subPlotColumns, subPlots(p));
        ax{p} = gca;
        yLimits(p, :) = ylim;
    end
    limits = [min(yLimits(:, 1)), max(yLimits(:, 2))];
    for p = 1:plots
        ylim(ax{p}, limits);
    end
end
