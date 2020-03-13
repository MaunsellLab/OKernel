function [IncludeDays] = makeIncludeDays()
%constructs a cell matrix of condition x animal where each cell contains
%all the dates to be included in okernel analysis

%the conditions are a step visual stimulus, a ramp visual stimulus, and a
%control condition

%load daily log
DailyLog = GetGoogleSpreadsheet('1-wv6DN9zflUfsD-HjmpE_YBHn4Omg1OjRR45bTaNJMg');

DateColumn = 1;
SubjectColumn = 2;
RampColumn = 7;
ControlFiberColumn = 18;
ControlVisStimColumn = 19;
VisibleKernelColumn = 23;

subjectList = unique(DailyLog(2:end,SubjectColumn));
subjectTotal = length(subjectList);
 
IncludeDays = cell(5,subjectTotal+1);
IncludeDays{1,1} = 'Mouse Number';
IncludeDays{2,1} = 'Step Stimulus'; StepRow = 2;
IncludeDays{3,1} = 'Ramp Stimulus'; RampRow = 3;
IncludeDays{4,1} = 'Control Opto Location';ControlOptoRow = 4;
IncludeDays{5,1} = 'Control Visual Stim Location'; ControlVisRow = 5;


for nowSubj = 1:subjectTotal
    IncludeDays{1,nowSubj+1} = subjectList{nowSubj};
    
    nowSubjNum = str2double(subjectList{nowSubj});
    
    for nowRow = 2:length(DailyLog(:,SubjectColumn))
        if str2double(DailyLog{nowRow,SubjectColumn}) == nowSubjNum
            if ~isempty(DailyLog{nowRow,VisibleKernelColumn})
                
                if isempty(DailyLog{nowRow,ControlFiberColumn}) && isempty(DailyLog{nowRow,ControlVisStimColumn})
                    ramp = str2double(DailyLog{nowRow,RampColumn});
                    
                    
                    if ramp == 0
                        IncludeDays{StepRow,nowSubj+1} = [IncludeDays{StepRow,nowSubj+1} ; DailyLog{nowRow,DateColumn}];
                    elseif ramp == 500
                        IncludeDays{RampRow,nowSubj+1} = [IncludeDays{RampRow,nowSubj+1} ; DailyLog{nowRow,DateColumn}];
                    else
                        disp(['not including ' DailyLog{nowRow,DateColumn} ' due to irregular ramp duration'])
                    end
                    
                    
                end
            end
            if ~isempty(DailyLog{nowRow,ControlFiberColumn})
                IncludeDays{ControlOptoRow,nowSubj+1} = [IncludeDays{ControlOptoRow,nowSubj+1} ; DailyLog{nowRow,DateColumn}];
            end
            if ~isempty(DailyLog{nowRow,ControlVisStimColumn})
                IncludeDays{ControlVisRow,nowSubj+1} = [IncludeDays{ControlVisRow,nowSubj+1} ; DailyLog{nowRow,DateColumn}];
            end
        end
    end % nowRow
    
    
    %%
    
end %nowSubj
 

end