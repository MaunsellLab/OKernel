function RTs = oneRTPerTrial(trials, indices)
%
% Return an array of RTs based on trials(:).reactTimeMS using logical array indices, ensuring that there is only 
% one RT for each trial. If more than one appears, use the first.

  RTs = [trials(indices).reactTimeMS];                       % get all the RTs
  if length(RTs) > sum(indices)                              % extra RTs?
    for t = 1:length(indices)
      if indices(t)                                          % selected trial?
        if length(trials(t).reactTimeMS) > 1                 % more than one RT?
          trials(t).reactTimeMS = trials(t).reactTimeMS(1);  % strip out extra RTs
        end
      end
    end
    RTs = [trials(indices).reactTimeMS];                     % remake the list
  end
end
