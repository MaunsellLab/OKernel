function grandAverageKernels
% Plot all the kernels for the grand average across selected sessions, once for steps and once for ramps

% 	dataDirName = '/Users/Shared/Data/OKernel/';
  [dataDirName, tableName] = whichData();
%   load([dataDirName ' Analysis/Mat Files/masterTable.mat'], 'T');
  load([dataDirName, tableName], 'T');
  limits = setLimits('All');
  rampLimits = [0, 500];
  for r = rampLimits
    limits.rampMS = r;
    U = selectUsingLimits(T, limits);
    if size(U, 1) == 0
      continue;
    end
    stimProfiles = getOptoProfiles(U);
    plotKernelPage(U, limits, stimProfiles);
    saveas(gcf, sprintf('%s Analysis/Figures/Kernels/Ramp %d %s.pdf', dataDirName, r, limits.animal{1}));
  end
end
