function [pvalue,varargout]=pvalueGet(y,input,method)  
if nargin==2
    method=1;% wenzhi's method; doing averaging across frames for each stimulus; method=2: do not do averaging;

end
if method==1
    angleNO=input.angleNO;
    trialNO=input.trialNO;
    framePerSti=input.framePerSti;
    y2=reshape(y,framePerSti,trialNO,angleNO);
    y3=squeeze(mean(y2,1));% first dimension is trialNO and the 2nd dimension is the angleNO
    [pvalue,~,~] = anova1(y3,[], 'off');
    if nargout>1
        y4=mean(y3);
        varargout{1}=y4;
    end
elseif method==2
    angleNO=input.angleNO;
    trialNO=input.trialNO;
    framePerSti=input.repeatNO;
    y2=reshape(y,framePerSti*trialNO,angleNO);
%     y3=squeeze(mean(y2,1));
    [pvalue,~,~] = anova1(y2,[], 'off');    
    
end