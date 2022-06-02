function [pvalue,y6,varargout]=pvalueGet_spatialTemproalMapping(y,input,method)  
if nargin==2
    method=1;% wenzhi's method; doing averaging across frames for each stimulus; method=2: do not do averaging;

end
if method==1
    angleNO=input.angleNO;
    trialNO=input.trialNO;
    framePerSti=input.framePerSti;
    %%
%     y2=reshape(y,framePerSti,trialNO,angleNO);
%     y3=squeeze(mean(y2,1));% first dimension is trialNO and the 2nd dimension is the angleNO
    y4=reshape(y,framePerSti,trialNO*angleNO);
    y4a=mean(y4(1:round(framePerSti/2),:),1);
    y4b=mean(y4(round(framePerSti/2+1):end,:),1);
    [y5,I]=max([y4a;y4b],[],1);
    y3=reshape(y5,trialNO,angleNO);
%     y3=reshape(y5,trialNO,angleNO);
%     y2=reshape(y,framePerSti,trialNO,angleNO);
%     y3=squeeze(mean(y2,1));% first dimension is trialNO and the 2nd dimension is the angleNO
    [pvalue,~,~] = anova1(y3,[], 'off');
    y6=mean(y3);
    if nargout>=3
        if mod(framePerSti,2)==1
            y4(end,:)=[];
            framePerSti=framePerSti-1;
        end            
            y_sd=y4(1:framePerSti/2,:);
            
            for ii=1:length(I)
                if I(ii)==2
                    y_sd(:,ii)=y4(framePerSti/2+1:end,ii);
                end
            end
            y_sd2=reshape(y_sd,framePerSti/2*trialNO,angleNO);
            ySD=std(y_sd2,0,1);
            ySE=ySD/sqrt(trialNO*framePerSti);
            varargout{1}=ySE;

        
    end
elseif method==2
    angleNO=input.angleNO;
    trialNO=input.trialNO;
    framePerSti=input.repeatNO;
    y2=reshape(y,framePerSti*trialNO,angleNO);
%     y3=squeeze(mean(y2,1));
    [pvalue,~,~] = anova1(y2,[], 'off');    
    y6=[];
end