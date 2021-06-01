function individualKernels
% Plot all the kernels for the grand average across selected sessions, once for steps and once for ramps

  rampMS = 0;
  animals = {'1462', '1463'};

  
  h = figure(2);
  set(h, 'Units', 'inches', 'Position', [25, 1.25, 8.5, 11]);
  clf;
  
  dataDirName = '/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/BehavData/';
  load([dataDirName ' Analysis/Mat Files/masterTable.mat'], 'T');
  limits = setLimits('All');
  limits.rampMS = rampMS;
  for a = 1:length(animals)
    limits.aniNum = a;
    limits.animal = animals{a};
    U = selectUsingLimits(T, limits);
    if height(U) == 0
      return;
    end
    stimProfiles = getOptoProfiles(U);
    plotKernelPage(U, limits, stimProfiles);
    saveas(gcf, sprintf('%s Analysis/Figures/Kernels/Ramp %d %s.pdf', dataDirName, limits.rampMS, limits.animal));
  end
  figure(2);
  sameAxisScaling('both', 4, 3, 1:length(animals));
  saveas(gcf, sprintf('%s Analysis/Figures/Kernels/Ramp %d Individuals.pdf', dataDirName, limits.rampMS));
end
