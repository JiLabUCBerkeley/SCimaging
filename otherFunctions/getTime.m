function [eyeT,timeRaw]=getTime(file)
    if nargin==0
        file='C:\serverData\2016.10.24 360961 nuclear\20161024_01*.png';
    end
fileName=getFileNames(file);
p=length(fileName);
file1=fileName(1).name;
eyeT=zeros(p,1);
timeRaw=zeros(p,6);
try
    dotNO=strfind(file1,'.');
    dotNO=dotNO(end);
    msFlag=strcmp(file1(dotNO-2:dotNO-1),'ms');% eyetracking files have to have ms at the end
    twoPflag=~(isempty(regexp(file1,'\d\d\d\d-\d\d-\d\d_','once')));
    stimulusFlag=~(isempty(regexp(file1,'_\d\d\d\d\d\d\d\d_\d\d\d\d\d\d','once')));
    if msFlag
        timeStamp=dotNO-3-7:dotNO-3;
        for ii=1:p
%             dotNO=strfind(file1,'.');
%             dotNO=dotNO(end);
%             timeStamp=dotNO-3-7:dotNO-3;
            eyeT(ii)=eval(fileName(ii).name(timeStamp));
        end
        eyeT=eyeT/1000;
    elseif twoPflag
         timeStamp=regexp(file1,'\d\d\d\d-\d\d-\d\d_')+(1:21)-1;
         formatIn='yyyy-mm-dd_HHMMSS_FFF';
         trash=file1(timeStamp);
         time0=datevec(trash,formatIn);
       for ii=1:p
           file1=fileName(ii).name;
           timeStamp=regexp(file1,'\d\d\d\d-\d\d-\d\d_')+(1:21)-1;
%            timeStamp=timeStamp2(1);
           timeII_string=file1(timeStamp);
           timeII=datevec(timeII_string,formatIn);
           eyeT(ii)=etime(timeII,time0);
           timeRaw(ii,:)=timeII;
       end
    elseif stimulusFlag
        timeStamp=regexp(file1,'_\d\d\d\d\d\d\d\d_\d\d\d\d\d\d')+(1:15);
        formatIn='yyyymmdd_HHMMSS';
        trash=file1(timeStamp);
        time0=datevec(trash,formatIn);
       for ii=1:p
           file1=fileName(ii).name;
           timeStamp=regexp(file1,'_\d\d\d\d\d\d\d\d_\d\d\d\d\d\d')+(1:15);
%            timeStamp=timeStamp2(1);
           timeII_string=file1(timeStamp);
           timeII=datevec(timeII_string,formatIn);
           eyeT(ii)=etime(timeII,time0);
           timeRaw(ii,:)=timeII;
       end        
    end    
catch me
end

function fileName=getFileNames(file)
starNO=strfind(file,'*');
if isempty(starNO)
    dotNO=strfind(file,'.');
    dotNO=dotNO(end);
    fileType=file(dotNO+1:end);
    fileName=dir([file(1:dotNO-1),'*.',fileType]);    
else
    fileName=dir(file);    
end

 