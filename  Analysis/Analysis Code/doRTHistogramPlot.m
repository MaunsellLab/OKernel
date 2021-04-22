 function doRTHistogramPlot(correctRTs, wrongRTs, missRTs, minRespTimeMS, maxRespTimeMS)
 % Plot the aggregated RT histogram
    subplot(4, 3, 2);
    cla;
    hold on;
    timeLimit = 750;
    RTBins = 100;
    while length(correctRTs) + length(missRTs) + length(wrongRTs) < RTBins * 10
      RTBins = RTBins / 2.0;
    end
    while true
      edges = linspace(-timeLimit, timeLimit, RTBins + 1);
      nCorrect = histc(correctRTs, edges); %#ok<*HISTC>
      nWrong = histc(wrongRTs, edges);
      nMiss = histc(missRTs, edges);
      if max([nWrong, nCorrect, nMiss] <= 50)         % re-bin on next plot?
          break;
      else
          RTBins = min([RTBins * 2, 100]);
      end
    end
    if sum(nCorrect) + sum(nWrong) + sum(nMiss) > 0
      binSize = edges(2) - edges(1);
      bH = bar(edges + binSize / 2, [nCorrect(:), nMiss(:), nWrong(:)], 'stacked');
      set(bH, 'barWidth', 1, 'lineStyle', 'none');
      set(bH(1), 'faceColor', [0 0 0.6]);
      set(bH(2), 'faceColor', [0.6 0 0]);
      set(bH(3), 'faceColor', [0.6 0 0]);
      yLimits = get(gca, 'YLim');                % vertical line at stimulus on
      plot([0 0], yLimits, 'k');
%       plot([minRespTimeMS, minRespTimeMS], yLimits, '--', 'Color', 0.5 * [0 1 0]);
%       plot([maxRespTimeMS, maxRespTimeMS], yLimits, '--', 'Color', 0.5 * [1 0 0]);
    end
    set(gca, 'XLim', [-timeLimit, timeLimit]);
    xlabel('Time Relative to Stimulus');
    title('Reaction Times');
end