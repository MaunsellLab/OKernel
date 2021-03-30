function displayDPrime
% Get a list of files to examine
  dataDirName = '/Users/Shared/Data/OKernel/';
	tableDataName = [dataDirName ' Analysis/Processed Files.mat'];
  limits.rampMS = 0;
  limits.criterion = 0;
  limits.animal = {'All'};
  limits.minTrials = 0;
  limits.minDec = -1;
  limits.oneDay = [];
  limits.minSessions = 0;
	[U, ~] = getSubset('normal', dataDirName, tableDataName, [], limits);
  U.Properties.VariableNames

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
%   U.pHit(find(U.pHit < 0)) = 0.5;
%   U.pFA(find(U.pFA > 0.9)) = 0.0;
  semilogy(U.dPrime, U.RTWindowMS, 'o');
  xlabel('d-prime');
  ylabel('Fit Response Window (ms)');
%   ylim([0, 0.5]);
  
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
%   ylim([0, 0.5]);
  
  animals = unique(U.animal);
  figure(4);
	set(gcf, 'units', 'inches', 'position', [27, 10.0, 7.5, 10]);    
  clf;
  pIndex = 1;
  for a = 1:length(animals)
    if strcmp(animals{a}, '905')
      fprintf('here');
    end
    thisAnimal = (U.animal == animals{a});
    if sum(thisAnimal) < 8
      continue;
    end
    subplot(4, 4, pIndex);
    histogram(U.noStimDPrime(U.animal == animals{a}) - U.stimDPrime(U.animal == animals{a}), 8);
    text(0.05, 0.95, {sprintf('d'' mean %.2f', ...
      nanmean(U.noStimDPrime(U.animal == animals{a}) - U.stimDPrime(U.animal == animals{a}))), ...
      sprintf('SD %.2f', std(U.dPrime(thisAnimal))), sprintf('n = %d', sum(thisAnimal))}, ...
      'units', 'normalized', 'horizontalAlignment', 'left', 'verticalAlignment', 'top');
    xlabel('delta d-prime');
    title(animals{a});
    pIndex = pIndex + 1;
  end
  sameXAxisScaling(4, 4, 1:pIndex - 1)
end