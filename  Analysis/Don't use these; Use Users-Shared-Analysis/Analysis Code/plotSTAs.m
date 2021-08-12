function plotSTAs

  load('STAsByUnit', 'unitSTA');
  figure(1);
  clf;
  for i = 1:size(unitSTA, 1)
    subplot(6, 7, i);
    plot(unitSTA(i,:));
    hold on;
    axis([0, size(unitSTA, 2), 0.3, 0.7]);
    plot([400, 400], [0.3, 0.7]);
  end
  
  numBS = 3;
  fprintf(sprintf('Running %d bootstrap samples', numBS));
  bootStrap = bootstrp(numBS, @mean, unitSTA);
  PCs = prctile(bootStrap, [5, 95]);
  
	figure(2);
  clf;
  plot(mean(unitSTA));
  hold on;
  CI05 = PCs(1,:);
  CI95 = PCs(2,:);
  x = 1:size(unitSTA, 2);
  x2 = [x, fliplr(x)];
  fillCI = [CI05, fliplr(CI95)];
  h = fill(x2, fillCI, 'b', 'lineStyle', '-', 'edgeColor', 'b', 'edgeAlpha', 0.5, 'faceAlpha', 0.10);
  set(h, 'faceAlpha', 0.10);
  ax = axis;
  plot([ax(1), ax(2)], [0.5, 0.5], 'k--');
  plot([400, 400], [ax(3), ax(4)], 'k--');
  title(sprintf('Weight by units (%d units, %d bootstrap samples)', size(unitSTA, 1), numBS));

end

