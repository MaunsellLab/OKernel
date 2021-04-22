function [dP, c] = dprime(pH, pFA, convertInf)
%
% Compute d-prime based on hits and false alarms 
%
% convert to Z scores and calculate d-prime and criterion
  if nargin < 3
    convertInf = false;
  end
  zH = -sqrt(2) .* erfcinv(2 * pH);
  zFA = -sqrt(2) .* erfcinv(2 * pFA);
  dP = zH - zFA ;
  c = -0.5 * (zH + zFA);
  if convertInf
    if abs(dP) == Inf
      dP = NaN;
    end
    if abs(c) == Inf
      c = NaN;
    end
  end
end
