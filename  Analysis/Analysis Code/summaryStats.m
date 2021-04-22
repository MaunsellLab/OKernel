function summaryStats
% For each animal, show the number of session and selected selection statistics

  [U, ~, ~] = getSessionTable('unfiltered');
  for r = [0, 500]
    fprintf('Ramp duration %d ms\n', r);
    V = U(U.rampMS == r, :);
    animals = unique(V.animal);
    [~, sortedIndices] = sort(str2double(animals));
    for a = sortedIndices'
      A = V(V.animal == animals{a}, :);
      fprintf(' %4s: %3d days; avg del-d'': %5.2f; mean power %5.3f\n', ...
        animals{a}, height(A), nanmean(A.noStimDPrime - A.stimDPrime), mean(A.meanPowerMW) ...
        );
    end
  end


end