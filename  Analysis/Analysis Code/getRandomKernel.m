function kernel = getRandomKernel(hits, misses, pulseDurMS, binsMS)
    % Get a random kernel. We initially tried using scrambled kernels from the experimental data, but that introduced
    % problems because hit and miss trials contain a lot of pre-selected structure, and because there were different
    % numbers of hits and missing with different effect sizes, scrambling didn't produce a expected value of 0.5.  Here
    % we simply create random kernels with the same structure as those used during the experiment.
    if hits + misses > 0
        hitKernel = getOneRandomdKernel(hits, pulseDurMS, binsMS);
        missKernel = getOneRandomdKernel(misses, pulseDurMS, binsMS);
        kernel = hitKernel - missKernel;
    else
        kernel = ones(1, binsMS) * 0.5;
    end
end

%%
function kernel = getOneRandomdKernel(numTrials, pulseDurMS, binsMS)
% Produce one binary white noise kernel based on randomly contructed stimuli
    if numTrials < 1 || pulseDurMS < 1 || binsMS < 1
    	kernel = ones(1, binsMS) * 0.5;
        return;
    end
    stimProfiles = zeros(numTrials, binsMS);
    for t = 1:numTrials
        profile = repelem(randi([0 1], 1, ceil(binsMS / pulseDurMS)), pulseDurMS);
        phaseBins = randi([0 (pulseDurMS - 1)]);                % get a random phase
        if phaseBins > 0
            prePend = repelem(randi([0 1]), phaseBins);
        else
            prePend = [];
        end
        stimProfiles(t, :) = [prePend profile(1:binsMS - phaseBins)];
    end
    kernel = mean(stimProfiles);
end