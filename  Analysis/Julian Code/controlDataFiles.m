function [dataFiles] = controlDataFiles(subjectnum,controlCondition)
%% only for plotted control data (does not set an inclusion criteria)
% inputs:
% subjectnum: subject number as a floating number. if more tan one, enclose
% in brackets
% controlCondition: 0 = pre-control; 1 = control; 2 = post-control;
% 3 = pre- and post-control combined; 4 is control with larger N to match
% combined pre and post control

% The plots should include:
%   1) Autocorrelogram of stimulus -- should be the triangle with 50 ms base (SEM for plots?)
%   2) Autocorrelogram of kernel -- should be broader than the stim auto-corr.
%   3) Cross-correlogram of stimulus-kernel (trace of N/P kernel?)

%datapath for development computer
datapath = '/Users/julian/Documents/MATLAB/OKernel/';
% outputpath = '/Users/julian/Documents/MATLAB/OKernel/';

%datapath for John's computer
% datapath = '/Users/Shared/Data/OKernel/';
outputpath = '/Users/Shared/Data/OKernel/';

%date ranges for each animal and condition

%1218 = Sophie
sophiePreCon  = {'2020-06-01','2020-06-05'}; %{'2020-05-31','2020-06-05'};
sophieCon = {'2020-06-12','2020-06-16'}; %{'2020-06-06','2020-06-16'}; %began controls 6/6
sophiePostCon = {'2020-06-21' , '2020-06-25'};%{'2020-06-17' , '2020-06-25'};

sophieSkip = {};
 
%1257 = Sufjan
sufjanPreCon = {'2020-05-25','2020-05-29'};
sufjanCon = {'2020-06-10','2020-06-14'}; %{'2020-06-10','2020-06-14'};%{'2020-05-30','2020-06-14'}; %began controls 5/30--6/7 and 6/9 not working?
sufjanPostCon = {'2020-06-21' , '2020-06-25'};%{'2020-06-15' , '2020-06-25'};% controlCondition: 0 = pre-control; 1 = control; 2 = post-control

%days to skip due to not enough trials
sufjanSkip = {'2020-06-08'};

%1220 = Caterina
caterinaPreCon = {'2020-06-17' '2020-06-21'};%{'2020-05-24' '2020-06-21'};
caterinaCon = {'2020-06-22' ,'2020-06-28'};
caterinaPostCon = {};

caterinaSkip = {'2020-06-16'};

%1150 = Joaquin
joaquinPreCon = {'2020-06-07' '2020-06-14'}; 
joaquinCon = {};
joaquinPostCon = {};

joaquinSkip = {'2020-06-04','2020-06-09', '2020-06-11','2020-06-12'};


if controlCondition == 3 || controlCondition == 4
    %control sessions must total 10 per animal b/c pre and post will each
    %be 5 per animal
    sophieCon = {'2020-06-07','2020-06-16'}; 
    sufjanCon = {'2020-06-04','2020-06-14'}; %skips over 06-08
    
end

if length(subjectnum) == 1
    nowsubject = num2str(subjectnum);
    
    if strcmp(nowsubject,'1218')
        
        skipDates = sophieSkip;
        
        if controlCondition == 0 || controlCondition == 3 %    pre-control
            daterange = sophiePreCon;
            if controlCondition == 3 
                daterange2 = sophiePostCon;
            end
        elseif controlCondition == 1 || controlCondition == 4 % control
            daterange = sophieCon;
        elseif controlCondition == 2
            daterange = sophiePostCon;
        else
            error('incorrect control condition: must be 0, 1, or 2')
        end
    elseif strcmp(nowsubject,'1257')
        
        skipDates = sufjanSkip;
        
        if controlCondition == 0 || controlCondition == 3  %     pre-control
            daterange = sufjanPreCon;
            if controlCondition == 3 
                daterange2 = sophiePostCon;
            end
        elseif controlCondition == 1 || controlCondition == 4 % control
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
        elseif controlCondition == 1 || controlCondition == 4 % control
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
        elseif controlCondition == 1 || controlCondition == 4 % control
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
    if controlCondition == 3 
        alldates2 = dateRange(daterange2{1},daterange2{2});
        alldates = [alldates alldates2];
    end
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
                [outputpath nowsubject '/MatFiles/' alldates{nowdate} '.mat']};
            
        end
    end
    
else %more than one subject number
    
    dataFiles = {};
    
    for nowsubj = 1:length(subjectnum)
        nowsubject = num2str(subjectnum(nowsubj));
        
        if strcmp(nowsubject,'1218')
            
            skipDates = sophieSkip;
            
            if controlCondition == 0 || controlCondition == 3  %    pre-control
                daterange = sophiePreCon;
                if controlCondition == 3 
                    daterange2 = sophiePostCon;
                end
            elseif controlCondition == 1 || controlCondition == 4 % control
                daterange = sophieCon;
            elseif controlCondition == 2
                daterange = sophiePostCon;
            else
                error('incorrect control condition: must be 0, 1, or 2')
            end
        elseif strcmp(nowsubject,'1257')
            
            skipDates = sufjanSkip;
            
            if controlCondition == 0 || controlCondition == 3 %     pre-control
                daterange = sufjanPreCon;
                if controlCondition == 3 
                    daterange2 = sufjanPostCon;
                end
            elseif controlCondition == 1 || controlCondition == 4 % control
                daterange = sufjanCon;
            elseif controlCondition == 2 % post-control
                daterange = sufjanPostCon;
            else
                error('incorrect control condition: must be 0, 1, or 2')
            end
        elseif strcmp(nowsubject,'1220')
            
            skipDates = caterinaSkip;
            
            if controlCondition == 0 || controlCondition == 3 %     pre-control
                daterange = caterinaPreCon;
            elseif controlCondition == 1 || controlCondition == 4% control
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
            elseif controlCondition == 1 || controlCondition == 4 % control
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
        
        if controlCondition == 3 && ~isempty(daterange2)
            alldates2 = dateRange(daterange2{1},daterange2{2});
            alldates = [alldates alldates2];
        end
        
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
                    [outputpath nowsubject '/MatFiles/' alldates{nowdate} '.mat']};
            end
        end
        daterange = [];
        daterange2 = [];
    end
end
subjectnum(isnan(subjectnum)) = [];
