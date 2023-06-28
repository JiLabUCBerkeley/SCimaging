function [result,output]=OSCriteria_moving_bars(handles,criteria)
if nargin==1  
    % Note: these criteria are overwritten by the criteria set in Step 4c 
    criteria.maxYThresh=10; % threshold to determine weather the cell is responsive; default value: 10 (10%)
    OS.annovaFlag=1;% option to choose whether include p value of annova test; 1: include the annova test; 0: does not include it;
    OS.fitFlag=0; % option to use fitted data or raw df/f to calculate OSI itself; 
    OS.goodFitFlag=[0.4, 0.6, 1]; % the first value is the threshold of SSE; the second value is the threshold of RSQ; the last one is option to use (1) or not to use (0) critetia of good fitting to designate the cell is OS;
    OS.gOSI_OSI=[0.25, 0.5, 1];% the first value is the threshold of gOSI; the second value is the threshold of OSI; the last one is to choose gOSI>0.25 (1) or OSI>0.5 (2)
    criteria.OS=OS;    
end

%% get criteria (input from step 4c)
% Default for considering an ROI to be OS is that max dff0 > 10%, pass
% anova, and good fit to gaussian
maxYThresh=criteria.maxYThresh; %max dff0 of ROI must be above threshold
OS=criteria.OS; %criteria for designating ROI's as OS (eg pass anova, good fit, etc)
annovaFlag=OS.annovaFlag; %if =1, fitted OSI rather than raw OSI is sued 
fitFlag=OS.fitFlag; %if =1, ROI must have good fit to gaussian
goodFitFlag=OS.goodFitFlag; %specifies SSE and RSQ required for good fit; if 3rd value=1 then ROI must have good fit
gOSI_OSI=OS.gOSI_OSI; %specifies whether ROI must have gOSI or OSI above some threshold (usually doesn't apply)

%% get values;
maxY=handles.maxY;
pvalue=handles.pvalue;
RSQ=handles.RSQ;
SSE=handles.SSE;
OSI_raw=handles.OSI_raw;
OSI_fitted=handles.OSI_fitted;
gOSI=handles.gOSI;

p=length(maxY); % #ROI

% Flag responsive ROIs (>10% dff0)
flag1=maxY>=maxYThresh; %logical vector with entry for each ROI; 1 means the ROI is responsive 
flag1b=flag1;
flag1c='OS_maxDf_f>10%'; %header for spreadsheet

%% whether to use the criterion that the pvalue should be smaller than 0.05 (annovaFlag(1))
flag2c='pvalue<0.05';
if annovaFlag(2)==1 %use anova criterion
    flag2b=pvalue<=annovaFlag(1);
    flag2=and(flag1,flag2b); %ROI is responsive and passes anova
    flag2c=[flag2c,', used'];
else
    flag2b=true(p,1);
    flag2=and(flag1,flag2b);
    flag2c=[flag2c,',not used'];
end
flag2c=['OS_', flag2c];

%% whether to use the criteria of good fitting 
flag3c='goodFitting';
if goodFitFlag(3)==1
    flag3b=SSE<goodFitFlag(1) & RSQ>goodFitFlag(2);
    flag3=and(flag2,flag3b); %ROI is responsive, passes anova, and is good fit
    flag3c=[flag3c,'SSE<',num2str(goodFitFlag(1),'%4.2f'),'RSQ>',num2str(goodFitFlag(2),'%4.2f')];
else
    flag3b=true(p,1);
    flag3=and(flag2,flag3b);
    flag3c=[flag3c,', not used'];
end

flag3c=['OS_', flag3c];
if fitFlag==0
    OSI=OSI_raw;
else
    OSI=OSI_fitted;
end


%% whether to use gOSI>0.25 (default value of gOSI_OSI(1)=0.25) as a criterion (when gOSI_OSI(3)=1) 
%or to use OSI>0.5 (default value of gOSI_OSI(2)=0.5), when gOSI_OSI(3)=2; typically not used as criteria 

if gOSI_OSI(3)==1 %use gOSI threshold
    flag4=gOSI>gOSI_OSI(1);
    flag4c=['gOSI>',num2str(gOSI_OSI(1),'%4.2f')];
elseif gOSI_OSI(3)==2 %use OSI threshold (which is typically 0; we just require good fit)
    flag4=OSI>gOSI_OSI(2);
    flag4c=['OSI>',num2str(gOSI_OSI(2),'%4.2f')];
    
        if fitFlag==0
            flag4c=['raw',flag4c];
        else
            flag4c=['fitted',flag4c];
        end   
end
flag4b=flag4;

% Make a "result" vector; logical vector specifying whether each ROI is OS (true) or not (false) 
result=and(flag3,flag4); %ROI is OS if it has good fit and passes gOSI/OSI threshold (usually 0 for OSI)


OS_flagAll=[result,flag1b,flag2b,flag3b,flag4b];
flagHead={'isOS?',flag1c,flag2c,flag3c,flag4c};
output.flagHead=flagHead;
output.OS_flagAll=OS_flagAll;
trash=0;


