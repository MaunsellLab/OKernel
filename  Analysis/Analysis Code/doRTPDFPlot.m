 function doRTPDFPlot(correctRTs, wrongRTs, missRTs, minRespTimeMS, maxRespTimeMS)

    subplot(4, 3, 3);
    cla;
    cdfplot([correctRTs wrongRTs missRTs]);
%     timeLimit = min(file.responseLimitMS, 5000);
    timeLimit = 1000;
    set(gca, 'XLim', [-500 timeLimit], 'YLim', [0 1]);
    hold on;
    yLimits = get(gca, 'YLim');
    plot([0 0], yLimits, 'k');
    plot(double(minRespTimeMS) * [1 1], yLimits, '--', 'Color', 0.5 * [0 1 0]);
    plot(double(maxRespTimeMS) * [1 1], yLimits, '--', 'Color', 0.5 * [1 0 0]);
%     plot(double(kernelRTMinMS) * [1 1], yLimits, ':', 'Color', 0.5 * [0 1 0]);
%     plot(double(kernelRTMaxMS) * [1 1], yLimits, ':', 'Color', 0.5 * [1 0 0]);
    xlabel('Time Relative to Stimulus');
    ylabel('');
    title('Cumulative Reaction Times');
 end