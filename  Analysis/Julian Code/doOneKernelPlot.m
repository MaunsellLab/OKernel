function doOneKernelPlot(subIndex, profile, startTimeMS, endTimeMS, plotTitle, yLabel, CI)
    subplot(4, 3, subIndex);
    cla;
  	bins = size(profile, 2);
    if CI > 0
        theMean = mean(profile(1:floor(bins/2)));               % take the mean as prestim period
        posCI = theMean + CI;
        negCI = theMean - CI;
        h = fill([0, bins, bins, 0], [posCI, posCI, negCI, negCI], [0.8, 0.8, 0.8]);
        set(h, 'linestyle', ':', 'facealpha', 0.25)
    end
    hold on;
    plot(profile, 'b');
    ax = gca;
    ax.XGrid = 'on';
    xlim(ax, [0, bins]);
    set(gca,'XTick', [0, -startTimeMS, -startTimeMS + 50 -startTimeMS + 100, bins]);
    set(gca, 'XTickLabel', {sprintf('%d', startTimeMS), '0', '50', '100', sprintf('%d', endTimeMS)});
    xlabel('Time Relative to Stimulus');
    ylabel(yLabel);
    title(plotTitle);
    hold off;
end
