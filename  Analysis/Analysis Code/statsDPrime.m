function statsDPrime
%
% Get the r^2 for the logistic fits to the RT distributions for steps and ramps and combined
%
  reportDPrime('All');
  reportDPrime('All Steps');
  reportDPrime('All Ramps');
end

function reportDPrime(subset)
  [U, ~] = getSessionTable(subset);
  fprintf('%s (n = %d)\n noStim: min d'' %.1f, max d'' %.1f\n', subset, height(U), min(U.noStimDPrime), max(U.noStimDPrime));
  fprintf('   stim: min d'' %.1f, max d'' %.1f\n', min(U.stimDPrime), max(U.stimDPrime));
  fprintf(' median base d'' %.1f (IQR %.1f - %.1f)\n', prctile(U.noStimDPrime, [50, 25, 75]));
  fprintf(' median delta d'' %.1f (IQR %.1f - %.1f)\n', prctile(U.noStimDPrime - U.stimDPrime, [50, 25, 75]));
end
