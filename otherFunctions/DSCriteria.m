function [result,output]=DSCriteria(handles,criteria)
if nargin==1
    
criteria.maxYThresh=10; % threshold to determine weather the cell is responsive; default value: 10 (10%)
DS.annovaFlag=1;% option to choDSe whether include p value of annova test; 1: include the annova test; 0: does not include it;
DS.fitFlag=0; % option to use fitted data or raw df/f to caluclate DSI itself; 
DS.goodFitFlag=[0.4, 0.6, 1]; % the first value is the threshold of SSE; the second value is the threshold of RSQ; the last one is option to use (1) or not to use (0) critetia of good fitting to designate the cell is DS;
DS.gDSI_DSI=[0.25, 0.5, 1];% the first value is the threshold of gDSI; the second value is the threshold of DSI; the last one is to choDSe gDSI>0.25 (1) or DSI>0.5 (2)
criteria.DS=DS;    
else

end

%% get criteria
maxYThresh=criteria.maxYThresh;
DS=criteria.DS;
annovaFlag=DS.annovaFlag;
fitFlag=DS.fitFlag;
goodFitFlag=DS.goodFitFlag;
gDSI_DSI=DS.gDSI_DSI;   

%% get values;
maxY=handles.maxY;
% dff0_avg=handles.dff0_avg;
pvalue=handles.pvalue;
RSQ=handles.RSQ;
SSE=handles.SSE;

% handles.DSI_raw=DSI_raw;
% handles.DSI_fitted=DSI_fitted;
DSI_raw=handles.DSI_raw;
DSI_fitted=handles.DSI_fitted;
gDSI=handles.gDSI;

p=length(maxY);
flag1=maxY>=maxYThresh;
flag1b=flag1;
flag1c='maxDf_f>10%';
flag1c=['DS_', flag1c];
%% whether to use the criterion that the pvalue should be smaller than 0.05 (annovaFlag(1))
flag2c='pvalue<0.05';
if annovaFlag(2)==1
    flag2b=pvalue<=annovaFlag(1);
    flag2=and(flag1,flag2b);
    flag2c=[flag2c,', used'];
else
    flag2b=true(p,1);
    flag2=and(flag1,flag2b);
    flag2c=[flag2c,',not used'];
end
flag2c=['DS_', flag2c];

%% whether to use the critetia of good fitting 
flag3c='goodFitting';
if goodFitFlag(3)==1
    flag3b=SSE<goodFitFlag(1) & RSQ>goodFitFlag(2);
    flag3=and(flag2,flag3b);
    flag3c=[flag3c,'SSE<',num2str(goodFitFlag(1),'%4.2f'),'RSQ>',num2str(goodFitFlag(2),'%4.2f')];
else
    flag3b=true(p,1);
    flag3=and(flag2,flag3b);
    flag3c=[flag3c,', not used'];
end

flag3c=['DS_', flag3c];
if fitFlag==0
    DSI=DSI_raw;
else
    DSI=DSI_fitted;
end


%% whether to use gDSI>0.25 (default value of gDSI_DSI(1)=0.25) as a standard (when gDSI_DSI(3)=1) 
%or to use DSI>0.5 (defaut value of gDSI_DSI(2)=0.5), when gDSI_DSI(3)=2; 

if gDSI_DSI(3)==1
    flag4=gDSI>gDSI_DSI(1);
    flag4c=['gDSI>',num2str(gDSI_DSI(1),'%4.2f')];
elseif gDSI_DSI(3)==2
    flag4=DSI>gDSI_DSI(2);
    flag4c=['DSI>',num2str(gDSI_DSI(2),'%4.2f')];
    
        if fitFlag==0
            flag4c=['raw',flag4c];
        else
            flag4c=['fitted',flag4c];
        end   
end
flag4b=flag4;

%% whether to use the condition that the cell is OS or not
if DS.OSFlag==0
    flag5b=true(p,1);
    flag5c='DS_ isOS? not used';
else
    flag5b=handles.OS.result;
     flag5c='DS_ isOS? used';
end

result=and(flag3,flag4);
result=and(result,flag5b);

DS_flagAll=[result,flag1b,flag2b,flag3b,flag4b,flag5b];
flagHead={'isDS?',flag1c,flag2c,flag3c,flag4c,flag5c};
output.flagHead=flagHead;
output.DS_flagAll=DS_flagAll;
trash=0;


