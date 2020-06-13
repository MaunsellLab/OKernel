%% outputs cell of dates between two input dates (yyyy-MM-dd)
function [alldates] = dateRange(firstdate,lastdate)

firstdatetime = datetime(firstdate,'InputFormat','yyyy-MM-dd');
lastdatetime = datetime(lastdate,'InputFormat','yyyy-MM-dd');

firstdatenum = datenum(firstdatetime);
lastdatenum = datenum(lastdatetime);

if lastdatenum < firstdatenum
    error('Date inputs must be chronological (first,last)')
end


alldates{1} = firstdate;
nowdate = firstdate;


datecounter = 1;
while strcmp(nowdate,lastdate) == 0
    datecounter = datecounter+1;
    nowdatetime = datetime(nowdate,'InputFormat','yyyy-MM-dd');
    nowdatetime = dateshift(nowdatetime, 'end','day');
    %     nowdate = dateshift(nowdate,'start','day','next');
    nowdate = datestr(nowdatetime,'yyyy-mm-dd');
    alldates{datecounter} = nowdate;
end


end

