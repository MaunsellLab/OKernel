function CI = stimCI(numTrials)
% 
% Return the 95% confidence interval for the binary optogenetic stimulus, normalized to a stimulus of 0 to 1.
% The variance of a summed binomial statistic is n*p*q, and SD (CI) is its square root. Dividing by n to bring
% the SD into the range of 0:1 gives sqrt(n*p*q)/n, or sqrt(p*q/n).  Because p=q=0.5 for our white noise stimuli,
% this becomes 0.5/sqrt(n). To convert the SD to 95% CI, we multiple by 1.96, assuming a normal distribution.

    CI = 1.96 * 0.5 / sqrt(numTrials);
end