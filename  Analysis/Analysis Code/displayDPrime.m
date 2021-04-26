function displayDPrime
% Get a list of files to examine

  rampMS = 0;
  animals = {'902', '1112', '1145', '905', '1223'};

%   rampMS = 500;
%   animals = {'902', '1112', '1150'};
  
	dataDirName = '/Users/Shared/Data/OKernel/';
  load([dataDirName ' Analysis/Mat Files/masterTable.mat'], 'T');
  limits = setLimits('All');
  limits.rampMS = rampMS;
	limits.animal = animals;
	U = selectUsingLimits(T, limits);
  if height(U) == 0
    return;
  end
  
  figure(3);
	set(gcf, 'units', 'inches', 'position', [27, 10.0, 7.5, 10]);    
  clf;
  subplot(3, 2, 1);
  histogram(U.dPrime, 25);
  text(0.05, 0.95, {sprintf('d'' mean %.2f', nanmean(U.dPrime(U.dPrime < 10000))), ...
    sprintf('SD %.2f', std(U.dPrime(U.dPrime < 10000)))}, 'units', 'normalized', ...
    'horizontalAlignment', 'left', 'verticalAlignment', 'top');
  xlabel('d-prime');
  
  subplot(3, 2, 2);
  plot(U.dPrime, U.pHit, 'o');
  hold on;
  plot(U.dPrime, U.pFA, 'o');
  xlabel('d-prime');
  ylabel('Probability');
  ylim([0, 1]);
  legend('p Hit', 'p False Alarm', 'location', 'northWest');

  subplot(3, 2, 3);
  plot(U.dPrime, U.c, 'o');
  xlabel('d-prime');
  ylabel('Criterion');
  a = axis;
  minAxis = min(a(1), a(3));
  maxAxis = max(a(2), a(4));
  axis([minAxis, maxAxis, minAxis, maxAxis]);
  hold on;
  
  subplot(3, 2, 4);
  semilogy(U.dPrime, U.RTWindowMS, 'o');
  xlabel('d-prime');
  ylabel('Fit Response Window (ms)');
  
  subplot(3, 2, 5);
  U.pHit(find(U.pHit < 0)) = 0.5;
  U.pFA(find(U.pFA > 0.9)) = 0.0;
  plot(U.pHit, U.pFA, 'o');
  xlabel('p Hit');
  ylabel('p False Alarm');
  ylim([0, 0.5]);
  
  subplot(3, 2, 6);
  U.pHit(find(U.pHit < 0)) = 0.5;
  U.pFA(find(U.pFA > 0.9)) = 0.0;
  plot(U.noStimDPrime - U.stimDPrime, U.kernelPeak, 'o');
  xlabel('opto change in d-prime');
  ylabel('kernel peak');
  
  animals = unique(U.animal);
  figure(4);
	set(gcf, 'units', 'inches', 'position', [27, 10.0, 7.5, 10]);    
  clf;
  pIndex = 1;
  for a = 1:length(animals)
    thisAnimal = (U.animal == animals{a});
    if sum(thisAnimal) < 8
      continue;
    end
    subplot(4, 4, pIndex);
    histogram(U.noStimDPrime(U.animal == animals{a}) - U.stimDPrime(U.animal == animals{a}), 8);
    text(0.05, 0.95, {sprintf('del-d'' mean %.2f', ...
      nanmean(U.noStimDPrime(U.animal == animals{a}) - U.stimDPrime(U.animal == animals{a}))), ...
      sprintf('SD %.2f', nanstd(U.dPrime(thisAnimal))), sprintf('n = %d', sum(thisAnimal))}, ...
      'units', 'normalized', 'horizontalAlignment', 'left', 'verticalAlignment', 'top');
    xlabel('delta d-prime');
    title(animals{a});
    pIndex = pIndex + 1;
  end
  sameXAxisScaling(4, 4, 1:pIndex - 1)
  
  figure(8);
	set(gcf, 'units', 'inches', 'position', [27, 10.0, 7.5, 10]);    
  clf;
  pIndex = 1;
  for a = 1:length(animals)
    thisAnimal = (U.animal == animals{a});
    if sum(thisAnimal) < 8
      continue;
    end
    U.noStimDPrime(thisAnimal & U.noStimDPrime == Inf) = NaN; 
    subplot(4, 4, pIndex);
    histogram(U.noStimDPrime(U.animal == animals{a}), 8);
    text(0.05, 0.95, {sprintf('no stim d'' mean %.2f', ...
      nanmean(U.noStimDPrime(U.animal == animals{a}))), ...
      sprintf('SD %.2f', nanstd(U.noStimDPrime(thisAnimal))), sprintf('n = %d', sum(thisAnimal))}, ...
      'units', 'normalized', 'horizontalAlignment', 'left', 'verticalAlignment', 'top');
    xlabel('d-prime');
    title(animals{a});
    pIndex = pIndex + 1;
  end
  sameXAxisScaling(4, 4, 1:pIndex - 1)

  goodAnimals = {'1223', '902', '905', '1112', '1145'};
  badAnimals = {'1218', '1220', '866', '1257', '1150', '1003', '1221', '844', '903'};

  figure(5);
  clf;
  peaksVDDPrime(1, goodAnimals, U, 'Good Animals');
  peaksVDDPrime(2, badAnimals, U, 'Bad Animals');
  sameAxisScaling('xy', 2, 1, 1:2);
  
  figure(6);
  clf;
  peaksVDPrime(1, goodAnimals, U, 'Good Animals');
  peaksVDPrime(2, badAnimals, U, 'Bad Animals');
  sameAxisScaling('xy', 2, 1, 1:2);
  
  figure(7);
  clf;
  str = cell(1, length(animals));
  for a = 1:length(animals)
    indices = (U.animal == animals{a});
    scatter(U.stimDPrime(indices), U.noStimDPrime(indices), 'filled');
    hold on;
    str{a} = sprintf('%s: %d', animals{a}, sum(indices));
  end
  text(0.05, 0.95, str, 'units', 'normalized', 'horizontalAlignment', 'left', 'verticalAlignment', 'top');
  xlabel('stim d-prime');
  ylabel('no stim d-prime');
  legend(animals, 'location', 'southeast');
  axis([-2, 4, -2, 4]);
  hold on;
  plot([-2, 4], [-2, 4], 'r-');
end

function peaksVDPrime(plotIndex, animals, U, theTitle)
  
  subplot(2, 1, plotIndex);
  str = cell(1, length(animals));
  for a = 1:length(animals)
    thisAnimal = (U.animal == animals{a});
    indices = thisAnimal & U.kernelPeak ~= 0;
    stimHitRate = double(U.stimCorrects(indices)) ./ double(U.stimCorrects(indices) + U.stimFails(indices));
    noStimHitRate = double(U.noStimCorrects(indices)) ./ double(U.noStimCorrects(indices) + U.noStimFails(indices));
    deltaDPrime = U.noStimDPrime(indices) - U.stimDPrime(indices);
    scatter(deltaDPrime, noStimHitRate - stimHitRate, 'filled');
    hold on;
    str{a} = sprintf('%s: %d', animals{a}, sum(indices & U.noStimDPrime - U.stimDPrime > 0.0));
  end
  text(0.8, 0.50, str, 'units', 'normalized', 'horizontalAlignment', 'left', 'verticalAlignment', 'top');
  xlabel('delta d-prime');
  ylabel('delta hit rate');
  title(theTitle);
  legend(animals);
end

function peaksVDDPrime(plotIndex, animals, U, theTitle)
  
  subplot(2, 1, plotIndex);
  str = cell(1, length(animals));
  for a = 1:length(animals)
    thisAnimal = (U.animal == animals{a});
    indices = thisAnimal & U.kernelPeak ~= 0;
    scatter(U.noStimDPrime(indices) - U.stimDPrime(indices), U.kernelPeak(indices), 'filled');
    hold on;
    str{a} = sprintf('%s: %d %.2f:%.2f', animals{a}, sum(indices), ...
      mean(U.noStimDPrime(indices) - U.stimDPrime(indices)), ...
      mean(U.kernelPeak(indices)));
  end
  text(0.8, 0.5, str, 'units', 'normalized', 'horizontalAlignment', 'left', 'verticalAlignment', 'top');
  xlabel('delta d-prime');
  ylabel('kernel peak (CIs)');
  title(theTitle);
  legend(animals);
end
