function [pvalue,y4]=pvalueGet_robustfit(y,input)  

    if input.bgMethod==2
        bcg=input.bcg;
    % if method==1
        angleNO=input.angleNO;
        trialNO=input.trialNO;
        framePerSti=input.framePerSti;
        y2=reshape(y,framePerSti,trialNO,angleNO);
        
%         brob=robustfit(bcg2,y2);
%         y2=y2-(brob(1)+brob(2)*bcg2);
        y3=squeeze(mean(y2,1));% first dimension is trialNO and the 2nd dimension is the angleNO
        y4=y3(:);bcg2=bcg(:);
        brob=robustfit(bcg2,y4);
        y4=y4-(brob(1)+brob(2)*bcg2);
        y3=reshape(y4,trialNO,angleNO);
        [pvalue,~,~] = anova1(y3,[], 'off');
        y4=mean(y3);        
    elseif input.bgMethod==1 || input.bgMethod==0
        angleNO=input.angleNO;
        trialNO=input.trialNO;
        framePerSti=input.framePerSti;
        y2=reshape(y,framePerSti,trialNO,angleNO);
        y3=squeeze(mean(y2,1));% first dimension is trialNO and the 2nd dimension is the angleNO
        [pvalue,~,~] = anova1(y3,[], 'off');
        y4=mean(y3);       
    elseif input.bgMethod==3
        pvalue=1; y4=nan;
    end
