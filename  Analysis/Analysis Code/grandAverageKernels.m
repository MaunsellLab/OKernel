function grandAverageKernels
% Plot all the kernels for the grand average across selected sessions, once for steps and once for ramps

  [dataDirName, tableName] = whichData();
  load([dataDirName, tableName], 'T');
  limits = setLimits('All');
  limits.rampMS = 0;
  U = selectUsingLimits(T, limits);
  if size(U, 1) == 0
    return;
  end
  stimProfiles = getOptoProfiles(U);
  plotKernelPage(U, limits, stimProfiles);
  saveas(gcf, sprintf('%s Analysis/Figures/Kernels/%s.pdf', dataDirName, limits.animal{1}));
end
