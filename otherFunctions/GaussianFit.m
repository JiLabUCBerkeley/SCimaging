function output=GaussianFit(test,order,output)
if nargin==2
    output.SEflag=1;
end
if isfield(output,'SEflag')
else
    output.SEflag=1;
end

[data3b,data3SEb,angleLeft2b]=avgRepeatNO(test(order(:,end),:),order(:,3),output);
%% Deal with negative values
minOptions=2; % 1 is wenzhi's way: move the whole curve to positive; 2 force negative values as zero; 0 means negative values will not be changed
output.minOptions=minOptions;
p=length(data3b);
if min(data3b)<0
    if minOptions==1
        minFlag=1;
        minY=min(data3b);
        data3b=data3b-min(data3b);

    elseif minOptions==2
        minFlag=1;
        ydataMinFlag=data3b(1:end-2)<0;
        if sum(ydataMinFlag)>=(p-2)/2
            output.minOptions=1;
            minY=min(data3b);
            data3b=data3b-min(data3b);
        else
            data3b(data3b<0)=0;  
            output.minOptions=2;
        end

    end

% else
%     minFlag=0;
%     output.minOptions=0; %negative values will remain negative
end
%%

[gDSI,pref_direction,exception]=gDSI_cal_Kath(data3b,length(data3b),output);
[gOSI,exception]=gOSI_cal_kath(data3b,length(data3b),output);

ii=1;
xx=1:size(data3b,1);xx=xx-1;xx=xx*360/size(data3b,1);
if nargin==2
    output.normFlag=1;
end
if isfield(output,'normFlag')
else
    output.normFlag=1;
end
% normFlag=input.normFlag;
% normFlag=1;
output=curveGaussianFitting_v2(xx,data3b(:,ii),data3SEb(:,ii),output);
output.x1=xx;
output.y1=data3b;
output.y1SE=data3SEb;
output.gDSI=gDSI;
output.gOSI=gOSI;
output.pref_direction=pref_direction;
output.maxY=max(data3b);
output.dff0_avg=data3b;