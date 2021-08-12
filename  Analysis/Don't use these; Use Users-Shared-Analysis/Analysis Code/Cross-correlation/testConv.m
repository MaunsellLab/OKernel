function testConv

  kernel = zeros(1, 800);
%   kernel(401:500) = (-1:-1:-100) / 100;
%   kernel(501:600) = (-100:-1) / 100;
  
%   kernel(301:400) = (-1:-1:-100) / 100;
  kernel(401:500) = (-100:-1) / 100;
  
  kernel = kernel + 0.05;                   % offset the kernel so the valid part of conv is visible

  figure(1);
  clf;
  sta = makeSTA(500, 350);
  c = conv(kernel, sta);
  doPlots(0, kernel, sta, c, 'conv(kernel, sta)', 'Conv.');              % simple convolution
  c = conv(kernel, sta, 'same');
  doPlots(1, kernel, sta, c, 'conv(kernel, sta, ''same'')', 'Conv.');    % force length to length(kernel)
  c = conv(sta, kernel, 'same');
  doPlots(4, kernel, sta, c, 'conv(sta, kernel, ''same'')', 'Conv.');    % forece length to length(sta)
  
  xTickLabels = {'-250', '-150', '-50', '50', '150', '250'};
  sta = makeSTA(500, 200);
  c = conv(kernel, sta, 'same');
  doPlots(2, kernel, sta, c, 'conv(kernel, sta, ''same'')', 'Conv.', xTickLabels);    % force length to length(kernel)   

  sta = makeSTA(500, 250);
  c = conv(kernel, sta, 'same');
  doPlots(3, kernel, sta, c, 'conv(kernel, sta, ''same'')', 'Conv.', xTickLabels);    % force length to length(kernel)   

  figure(2);
  clf;
  sta = makeSTA(500, 350);

  c = xcorr(kernel, sta);
  doPlots(0, kernel, sta, c, 'xcorr(kernel, sta)', 'XCorr.');              % simple convolution
  c = xcorr(kernel, sta, length(kernel) / 2 - 1);
  doPlots(1, kernel, sta, c, 'xcorr(kernel, sta, maxlag)', 'XCorr.');    % force length to length(kernel)
  c = xcorr(sta, kernel);
  doPlots(4, kernel, sta, c, 'xcorr(sta, kernel)', 'XCorr.');    % forece length to length(sta)
  
  xTickLabels = {'-400', '-300', '-200', '-100', '0', '100'};
  sta = makeSTA(500, 400);
%   c = xcorr(kernel, sta, length(kernel) / 2 - 1);
  c = xcorr(kernel, sta);
  doPlots(2, kernel, sta, c, 'xcorr(kernel, sta)', 'XCorr.', xTickLabels);    % force length to length(kernel)   

  sta = makeSTA(500, 300);
%   c = xcorr(kernel, sta, length(kernel) / 2 - 1);
  c = xcorr(kernel, sta);
  doPlots(3, kernel, sta, c, 'xcorr(kernel, sta)', 'XCorr.', xTickLabels);    % force length to length(kernel)   

end

function doPlots(row, kernel, sta, c, subtitle, funcStr, xTickLabels)
  numRows = 5;
  numCols = 3;
  
  subplot(numRows, numCols, row * numCols + 1);
  plot(kernel);
  hold on;
  a = axis;
  plot([400, 400], [a(3), a(4)], 'k:');
  set(gca,'XTick', 0:200:800);
  set(gca, 'XTickLabel', {'-400', '-200', '0', '200', '400'});
  [~, minIndex] = min(kernel);
  title(sprintf('Kernel (Peak at %d)\n(100 ms post-stim)', minIndex));
  
  subplot(numRows, numCols, row * numCols + 2)
  cla;
  plot(sta);
  hold on;
  a = axis;
  a(2) = length(sta);
  axis(a);
  set(gca,'XTick', 0:100:500);
  if nargin > 6
    set(gca, 'XTickLabel', xTickLabels);
    plot([400, 400], [a(3), a(4)], 'k:');
  else
    set(gca, 'XTickLabel', {'-400', '-300', '-200', '-100', '0', '100'});
    plot([400, 400], [a(3), a(4)], 'k:');
  end
  [~, minIndex] = min(sta);
  title(sprintf('STA (Peak at %d)', minIndex));
  
  subplot(numRows, numCols, row * numCols + 3)
  plot(1:length(c), c);
  axis([0, length(c), -inf, inf]);
  [~, maxIndex] = max(c);
  title({sprintf('%s (Peak at %d)', funcStr, maxIndex), subtitle});
	set(gca,'XTick', [1, length(c)]);
  set(gca, 'XTickLabel', {'1', sprintf('%d', length(c))});
end

function	sta = makeSTA(len, deltaIndex)
  sta = zeros(1, len);
  sta(deltaIndex) = -1.0;
end
