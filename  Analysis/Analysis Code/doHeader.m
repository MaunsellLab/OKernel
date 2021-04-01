function doHeader(U, limits)
%
% Display the header information for an OKernel analysis page
%
  headerText = cell(1, 1);
  if limits.rampMS == 0
    headerText{1} = sprintf('Visual Contrast Step');
  else
    headerText{1} = sprintf('Visual Contrast Ramp %d ms', limits.rampMS);
  end
  if strcmp(limits.animal, 'All')
    headerText{end + 1} = sprintf('%d sessions from %d animals', size(U, 1), length(unique(U.animal)));
    animalList = [];
    animals = unique(U.animal);
    for a = 1:length(animals)
      animalList = strcat(animalList, sprintf('  %s (%d),', animals{a}, sum(U.animal == animals{a})));
      if ~mod(a, 4)
        headerText{end + 1} = animalList; %#ok<AGROW>
        animalList = [];
      end
    end
    if mod(a, 4)
      headerText{end + 1} = animalList; 
    end
  else
    headerText{end + 1} = sprintf('%d sessions from Animal %s', size(U, 1), limits.animal);
  end
  headerText{end + 1} = sprintf('%d trial minimum for hits/misses', limits.minTrials);
  if limits.minDec == -1
    headerText{end + 1} = 'No required delta-d'' with opto';
  else
    headerText{end + 1} = sprintf('Delta d'' >=%.2f', limits.minDec);
  end
  headerText{end + 1} = sprintf('Avg(SD) d'' noStim: %.2f (%.2f)', mean(U.noStimDPrime), std(U.noStimDPrime));
  headerText{end + 1} = sprintf('Avg(SD) d''   stim: %.2f (%.2f)', mean(U.stimDPrime), std(U.stimDPrime));
  headerText{end + 1} = sprintf('Avg(SD) delta d'': %.2f (%.2f)', mean(U.noStimDPrime - U.stimDPrime), ...
    std(U.noStimDPrime - U.stimDPrime));
  
  axisHandle = subplot(4, 3, 1);						% default axes are 0 to 1
  set(axisHandle, 'visible', 'off');
  set(axisHandle, 'outerPosition', [0.02 0.75, 0.25, 0.2]);
  text(0.00, 1.25, 'OKernel', 'FontWeight', 'bold', 'FontSize', 16);
  text(0.00, 1.10, headerText, 'VerticalAlignment', 'top');
end
