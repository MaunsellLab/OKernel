function result = bfd(yData, pulseWidth)
% Brute force deconvolution modeling
%
  if nargin == 0
    testBFD();
    return;
  end

  numSamples = length(yData);
%   boundsFactor = 0.045;
%   lb = yData - boundsFactor;
%   ub = yData + boundsFactor;
  result  = lsqcurvefit(@myFun, yData, numSamples, yData);
  
    % Nested function (to get access to pulse width)
    function result = myFun(x, ~)
      f1 = zeros(1, numSamples);
      f1((1:pulseWidth) + numSamples / 2) = 1.0;
      result = conv(x, f1, 'same');
%       result = conv(x, f1, 'same') + 100 * x;
    end
end

function testBFD

  numSamples = 800;
  peak = 1.0;
  peakWidth = 10;
  pulseWidth = 25;
%   peakInc = peak / peakWidth;
%   noise = 0.1;
  numCols = 2;
  numRows = 2;
  
  figure(1);
	set(gcf, 'units', 'inches', 'position', [27, 10.0, 7.5, 10]);    
  clf;

% 	x = zeros(1, numSamples);
%   x(numSamples / 2 - peakWidth:numSamples / 2) = 0.0:peakInc:peak;
%   x(numSamples / 2 + 1:numSamples / 2 + peakWidth) = peak - peakInc:-peakInc:0.0;
%   x = x + rand(1, numSamples) * noise - noise / 2.0;
%   subplot(numCols, numRows, 1);
%   plot(x);
%   title('Raw Signal, Pre-convolution');
  
  f = zeros(1, numSamples);
  f([1:pulseWidth] + numSamples / 2) = 1.0;
  subplot(numCols, numRows, 2);
  plot(f);
  title('Pulse to Deconvolve');
  
%   c = conv(x, f, 'same');
  load('/Users/maunsell/Desktop/Step Kernel.mat', 'CIs'); 
  c = CIs(2, :);
  meanC = mean(c);
  c = c - meanC;
  subplot(numCols, numRows, 1);
  plot(c);
  title('Neuronal-Behavior Kernel');
  
  result = bfd(c, pulseWidth);
  subplot(numCols, numRows, 3);
  plot(result);
  title('Deconvolved Kernel');
  
  convResult = conv(result, f, 'same');
  subplot(numCols, numRows, 4);
  plot(convResult);
	title('Deconvolved Kernel Convolved with Pulse');

  
  sameAxisScaling('x', numCols, numRows, [1, 2, 4]);

end
