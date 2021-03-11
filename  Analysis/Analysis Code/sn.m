function sn

length = 100;
for t = 1:3
  trials = 10^(t+1);
  a = rand(trials, length);
  noise = std(mean(a));
  a(1, length/2) = 50;
  meanA = mean(a);
  peakA = max(meanA);
  signal = peakA - mean(meanA);

  subplot(3, 1, t);
  plot(meanA);
  ylim([0, 1]);
  text(0.05, 0.95, sprintf('%d trials\nS: %.3f\nN: %.3f\nSN: %.1f', trials, signal, noise, signal / noise), ...
    'units', 'normalized', 'verticalAlignment', 'top');
  if t == 1
    title('Fixed signal trials, varying noise trials');
  end
end