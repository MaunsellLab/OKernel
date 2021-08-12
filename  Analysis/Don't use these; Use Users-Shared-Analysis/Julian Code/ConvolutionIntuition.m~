

%% autocorrelated deconvolution with a delta function kernel
% box = zeros(1,100);
% box(51:75) = 1;
box = ones(1,25);

kdelta = zeros(1,1000);
kdelta(550) = 1;

kconvolution = conv(kdelta,box,'same');

figure; hold on;
plot(kconvolution);
plot(kdelta,'r')
legend('signal convolved with 25 ms boxcar','original signal')
title('Delta convolved with 25 ms box');
hold off;

kdeconvolution = deconv(kconvolution,box);

figure; hold on;
plot(kdeconvolution);
plot(kdelta,'r')
legend('deconvolved signal using 25 ms boxcar','original signal')
title('Convolved function then deconvolved with 25 ms box');
hold off;


%% temporal precision analysis for for STA x kernel convolution

boxcar = ones(1,25);
x1 = -199:1:200;

kernelreal = -rescale(gaussmf(x1,[5 60]));
[toss,kernelreal_peakidx] = min(kernelreal);

stareal = rescale(gaussmf(x1,[10 -50]) - gaussmf(x1,[25 -25]),-1,.5);
[toss,stareal_peakidx] = min(stareal);

staconvbox = rescale(conv(stareal,box,'same'),min(stareal),max(stareal));
kernelconvbox = rescale(conv(kernelreal,box,'same'),min(kernelreal),max(kernelreal));

staxkernelreal = conv(kernelreal,stareal,'same');
staxkernelbox = conv(kernelconvbox,staconvbox,'same');
[toss,staxkernelreal_peakidx] = max(staxkernelreal);

staxcorr = xcorr(kernelreal,stareal,200);
[toss,staxcorr_peakidx] = max(staxcorr);

figure; hold on;

subplot(3,1,1); hold on;
title(' "real" kernel conv w/ boxcar')
plot(x1,kernelreal,'k-')
plot(x1,kernelconvbox,'b-')
line([0 0],[-1 1],'LineStyle',':')
xlim([x1(1) x1(end)])

subplot(3,1,2); hold on;
title(' "real" sta conv w/ boxcar')
plot(x1,stareal,'k-')
plot(x1,staconvbox,'b-')
line([0 0],[-1 1],'LineStyle',':')
xlim([x1(1) x1(end)])

subplot(3,1,3); hold on;
title(' "real" sta conv w/ "real krnel" vs. conv w/ boxcar')
plot(x1,staxkernelreal,'k-')
plot(x1,staxkernelbox,'b-')
line([0 0],[-1 1],'LineStyle',':')
xlim([x1(1) x1(end)])



%% convolution vs. cross-correlation for STA x kernel

boxcar = ones(1,25);
x1 = -199:1:200;
x2 = -199:1:200;

kernelreal = -rescale(gaussmf(x1,[10 60]));
[toss,kernelreal_peakidx] = min(kernelreal);

stareal = rescale(gaussmf(x2,[10 -50]) - gaussmf(x2,[25 -25]),-1,.5);
[toss,stareal_peakidx] = min(stareal);

staconvbox = rescale(conv(stareal,box,'same'),min(stareal),max(stareal));
kernelconvbox = rescale(conv(kernelreal,box,'same'),min(kernelreal),max(kernelreal));

staxkernelreal = conv(kernelreal,stareal,'same');
staxkernelbox = conv(kernelconvbox,staconvbox,'same');
[toss,staxkernelreal_peakidx] = max(staxkernelreal);

staxcorr = xcorr(kernelreal,stareal,200);
[toss,staxcorr_peakidx] = max(staxcorr);

figure; hold on;

subplot(1,3,1); hold on;
title('Kernel')
plot(x1,kernelreal)
text(-180,.8,['kernel peak = ' num2str(x1(kernelreal_peakidx)) ' ms'])
line([0 0],[-1 1],'LineStyle',':')
xlim([x1(1) x1(end)])

subplot(1,3,2); hold on;
title('STA')
plot(x2,stareal)
text(-180,.8,['STA peak = ' num2str(x1(stareal_peakidx)) ' ms'])
line([0 0],[-1 1],'LineStyle',':')
xlim([x1(1) x1(end)])

subplot(1,3,3); hold on;
title('x-corr (red) vs conv (blue)')
plot(x1,staxcorr(1:end-1),'r-')
plot(x1,staxkernelreal,'b-')
text(-190,15,['X-corr peak = ' num2str(x1(staxcorr_peakidx)) ' ms'])
text(-190,12,['Conv peak = ' num2str(x1(staxkernelreal_peakidx)) ' ms'])
line([0 0],[-10 20],'LineStyle',':')
xlim([x1(1) x1(end)])