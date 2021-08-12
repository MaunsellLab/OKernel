function [plotStartMS, plotEndMS, plotRTStartMS] = plotLimits(onLine)

if nargin > 0 && onLine
  plotStartMS = -200;
  plotEndMS = 200;
  plotRTStartMS = plotStartMS - 200;
else
  plotStartMS = -400;
  plotEndMS = 400;
  plotRTStartMS = plotStartMS - 200;
end
