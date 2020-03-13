% Animates kernel construction from a single day's data


SyncData

nowDate = '2019-10-10';


subj = '902';
datapath =  '/Users/julian/Documents/MATLAB/OKernel/';
analysispath = '/Users/julian/Documents/MATLAB/OKernel/Analysis/Julian Code';
datafolder = [datapath subj '/MatFiles'];
    
cd(datafolder)

nowfile = [nowDate '.mat'];

load(nowfile)

% OKMatlab_191030(dParams, file, trials);
meanPowers = [trials(:).meanPowerMW];
stimIndices = meanPowers > 0;
trialEnds = [trials(:).trialEnd];
hitIndices = stimIndices & trialEnds == 0;
missIndices = stimIndices & trialEnds == 2;
numStimtrials = sum([hitIndices missIndices]);
    
for tr = 1:length(trials)
    nowtrials = trials(1:tr);
    if trials(tr).meanPowerMW == 0
        continue
    end
plotDayKernel(~, file, nowtrials)
end
