%{

Analysis code overview:

statsRTFit produces statistics on the analysis of response windows.  These involve fitting a logistic function
to the RT cumulative distribution and assigning a response window to the 2.5%-97.5% portion of the sigmoid.
The values reported are the median and IQR for the response window width for steps, ramps and all, together
with a ranksum test between the widths for steps and ranks.

testDPrime was written to test the performance of the d-prime calculation.  It also produces a separate figure for
each file that shows the RT distribution, the cumulative RT distribution, and the fitted response window.

%}
