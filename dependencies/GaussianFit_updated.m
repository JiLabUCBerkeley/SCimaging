function output=GaussianFit_updated(data,order,output)
if nargin==2
    output.SEflag=1;
end
if isfield(output,'SEflag')
else
    output.SEflag=1;
end

% data(data<0)=0; %set negative values to zero (assuming they're noise)
% output.minOptions=1;

[data3b,data3SEb,angleLeft2b]=avgRepeatNO_updated(data(order(:,end),:),order(:,3),output);
    %test(order(:,end),:) sorts the dff0 values by the 4th col of array
        %"order," which is the frame indices of the diff stimuli in their original
        %sequence
    %order(:,3) is the 3rd col of "order," which is the sorted stimulus
        %angles
    %so the inputs to the avgRepeatNO fxn are dff0 values sorted by the
        %stim angle, eg the values for 0 deg come first, then values for 30
        %deg, etc
    %output data3b is the average dff0 values for each stim direction,
        %avg'd across stimulus repeats and frames, so you have one value
        %per stimulus
    %output data3SEb is standard error of mean
    %output angleLeft2b is a vector with values for all stim angles (e.g.
        %0-330)
    
%% Deal with negative values
minOptions=output.negativeMethod; %1 to set negative values to zero; 2 to shift all values upward; 0 means negative values will not be changed

output.minOptions=minOptions;
p=length(data3b);
if min(data3b)<0
    if minOptions==1
        data3b(data3b<0)=eps; %set the negative values to ~zero (this assumes the negative values are noise)
        output.minOptions=1;
    elseif minOptions==2
        data3b=data3b-min(data3b); %increase all the values
        output.minOptions=2;
    end
else  %no negative values in data; don't change anything
    output.minOptions=0; 
end
%%
[gDSI,pref_direction,exception]=gDSI_cal_updated(data3b,length(data3b),output);
[gOSI,exception]=gOSI_cal_updated(data3b,length(data3b),output);

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
output=curveGaussianFitting_v2_updated(xx,data3b(:,ii),data3SEb(:,ii),output);
output.x1=xx;
output.y1=data3b; %this is the same as dff0_avg
output.y1SE=data3SEb;
output.gDSI=gDSI;
output.gOSI=gOSI;
output.pref_direction=pref_direction;
output.maxY=max(data3b);
output.dff0_avg=data3b; %dff0 is avg'd across stimulus repeats and frames