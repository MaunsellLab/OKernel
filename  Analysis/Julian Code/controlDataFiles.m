function [dataFiles] = controlDataFiles(subjectnum,controlCondition)
%% only for plotted control data (does not set an inclusion criteria)
% inputs:
% subjectnum: subject number as a floating number. if more tan one, enclose
% in brackets
% controlCondition: 0 = pre-control; 1 = control; 2 = post-control

% The plots should include:
%   1) Autocorrelogram of stimulus -- should be the triangle with 50 ms base (SEM for plots?)
%   2) Autocorrelogram of kernel -- should be broader than the stim auto-corr.
%   3) Cross-correlogram of stimulus-kernel (trace of N/P kernel?)

%datapath for development computer
% datapath = '/Users/julian/Documents/MATLAB/OKernel/';

%datapath for John's computer
datapath = '/Users/Shared/Data/OKernel/';


%date ranges for each animal and condition

%1218 = Sophie
sophiePreCon  = {'2020-05-31','2020-06-05'};
sophieCon = {'2020-06-06','2020-06-13'}; %began controls 6/6
sophiePostCon = {};

sophieSkip = {};
 
%1257 = Sufjan
sufjanPreCon = {'2020-05-25','2020-05-31'};
sufjanCon = {'2020-05-30','2020-06-13'}; %began controls 5/30--6/7 and 6/9 not working?
sufjanPostCon = {};% controlCondition: 0 = pre-control; 1 = control; 2 = post-control

%days to skip due to not enough trials
sufjanSkip = {'2020-06-08'};

%1220 = Caterina
caterinaPreCon = {'2020-05-24' '2020-06-13'};%{'2020-05-24' '2020-06-11'};
caterinaCon = {};
caterinaPostCon = {};

caterinaSkip = {};

%1150 = Joaquin
joaquinPreCon = {'2020-06-07' '2020-06-13'}; 
joaquinCon = {};
joaquinPostCon = {};

joaquinSkip = {'2020-06-04','2020-06-09', '2020-06-11','2020-06-12'};

if length(subjectnum) == 1
    nowsubject = num2str(subjectnum);
    
    if strcmp(nowsubject,'1218')
        
        skipDates = sophieSkip;
        
        if controlCondition == 0 %    pre-control
            daterange = sophiePreCon;
        elseif controlCondition == 1 % control
            daterange = sophieCon;
        elseif controlCondition == 2
            daterange = sophiePostCon;
        else
            error('incorrect control condition: must be 0, 1, or 2')
        end
    elseif strcmp(nowsubject,'1257')
        
        skipDates = sufjanSkip;
        
        if controlCondition == 0 %     pre-control
            daterange = sufjanPreCon;
        elseif controlCondition == 1 % control
            daterange = sufjanCon;
        elseif controlCondition == 2 % post-control
            daterange = sufjanPostCon;
        else
            error('incorrect control condition: must be 0, 1, or 2')
        end
    elseif strcmp(nowsubject,'1220')
        
        skipDates = caterinaSkip;
        
        if controlCondition == 0 %     pre-control
            daterange = caterinaPreCon;
        elseif controlCondition == 1 % control
            daterange = caterinaCon;
        elseif controlCondition == 2 % post-control
            daterange = caterinaPostCon;
        else
            error('incorrect control condition: must be 0, 1, or 2')
        end
    elseif strcmp(nowsubject,'1150')
        
        skipDates = joaquinSkip;
        
        if controlCondition == 0 %     pre-control
            daterange = joaquinPreCon;
        elseif controlCondition == 1 % control
            daterange = joaquinCon;
        elseif controlCondition == 2 % post-control
            daterange = joaquinPostCon;
        else
            error('incorrect control condition: must be 0, 1, or 2')
        end
    end
    
    if strcmp(daterange{2},'today')
        daterange{2} = datestr(datetime('today'),'yyyy-mm-dd');
    end
    
    alldates = dateRange(daterange{1},daterange{2}); %must be format yyyy-mm-dd
    
    for i = 1:length(skipDates)
        for j = 1:length(alldates)
            if strcmp(skipDates{i},alldates{j})
                alldates{j} = [];
            end
        end
    end
    

    
    dataFiles = {};
    
    
    for nowdate = 1:length(alldates)
        if isfile([datapath nowsubject '/MatFiles/' alldates{nowdate} '.mat'])
            
            dataFiles = {dataFiles{:},...
                [datapath nowsubject '/MatFiles/' alldates{nowdate} '.mat']};
            
        end
    end
    
else %more than one subject number
    
    dataFiles = {};
    
    for nowsubj = 1:length(subjectnum)
        nowsubject = num2str(subjectnum(nowsubj));
        
        if strcmp(nowsubject,'1218')
            
            skipDates = sophieSkip;
            
            if controlCondition == 0 %    pre-control
                daterange = sophiePreCon;
            elseif controlCondition == 1 % control
                daterange = sophieCon;
            elseif controlCondition == 2
                daterange = sophiePostCon;
            else
                error('incorrect control condition: must be 0, 1, or 2')
            end
        elseif strcmp(nowsubject,'1257')
            
            skipDates = sufjanSkip;
            
            if controlCondition == 0 %     pre-control
                daterange = sufjanPreCon;
            elseif controlCondition == 1 % control
                daterange = sufjanCon;
            elseif controlCondition == 2 % post-control
                daterange = sufjanPostCon;
            else
                error('incorrect control condition: must be 0, 1, or 2')
            end
        elseif strcmp(nowsubject,'1220')
            
            skipDates = caterinaSkip;
            
            if controlCondition == 0 %     pre-control
                daterange = caterinaPreCon;
            elseif controlCondition == 1 % control
                daterange = caterinaCon;
            elseif controlCondition == 2 % post-control
                daterange = caterinaPostCon;
            else
                error('incorrect control condition: must be 0, 1, or 2')
            end
        elseif strcmp(nowsubject,'1150')
            
            skipDates = joaquinSkip;
            
            if controlCondition == 0 %     pre-control
                daterange = joaquinPreCon;
            elseif controlCondition == 1 % control
                daterange = joaquinCon;
            elseif controlCondition == 2 % post-control
                daterange = joaquinPostCon;
            else
                error('incorrect control condition: must be 0, 1, or 2')
            end
        end
        
        if isempty(daterange)
            subjectnum(nowsubj) = NaN;
            
            continue
        end
        
        if strcmp(daterange{2},'today')
            daterange{2} = datestr(datetime('today'),'yyyy-mm-dd');
        end
        
        alldates = dateRange(daterange{1},daterange{2}); %must be format yyyy-mm-dd
        
        for i = 1:length(skipDates)
            for j = 1:length(alldates)
                if strcmp(skipDates{i},alldates{j})
                    alldates{j} = [];
                end
            end
        end
    
        for nowdate = 1:length(alldates)
            if isfile([datapath nowsubject '/MatFiles/' alldates{nowdate} '.mat'])
            
                dataFiles = {dataFiles{:},...
                    [datapath nowsubject '/MatFiles/' alldates{nowdate} '.mat']};
            end
        end
    end
end