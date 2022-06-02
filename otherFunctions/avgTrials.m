function [data3,data3Std,angleLeft2]=avgTrials(data2,angleLeft,repeatNO,angleNO,trialNO,input)
qt=[trialNO,angleNO];
[m,n,o]=size(data2);

if o>1
    o2=floor(o/trialNO);
    angleLeft2=zeros(o2,1);
    [Y,I]=sort(angleLeft);
    t=I;
    ii=1;
    t2=reshape(t,repeatNO,qt(ii,1),qt(ii,2));
    jjTmp=0;
    t3=cell(o2,1);
    for jj1=1:qt(ii,2)% how many angles
        for jj2=1:repeatNO
            jjTmp=jjTmp+1;
            chosenI=t2(jj2,:,jj1);
            t3{jjTmp}=squeeze(chosenI);

        end
    end

    data3=zeros(m,n,o2,'uint16');
    data3Std=data3;        
     for jj=1:o2
        imgRaw=double(data2(:,:,t3{jj}));
    %         [imgRaw,shiftAll,max2]=imageRegistration(imgRaw,1);
    %         [imgRaw,shiftAll,max2]=imageRegistration_previousRef(imgRaw,1,1);
    %         [imgAll2,shiftAll,max2]=imageRegistration_previousRef(imgAll,upsampling,firstRef)
        img=nanmean(imgRaw,3);
        img2=uint16(img);
        data3(:,:,jj)=img2;
        data3Std(:,:,jj)=uint16(nanstd(imgRaw,0,3));
        angleLeft2(jj)=angleLeft(t3{jj}(1));



    end    
else
    o=m;
    o2=floor(o/trialNO);
    angleLeft2=zeros(o2,1);
    [Y,I]=sort(angleLeft);
    t=I;
    ii=1;
    t2=reshape(t,repeatNO,qt(ii,1),qt(ii,2));
    jjTmp=0;
    t3=cell(o2,1);
    for jj1=1:qt(ii,2)% how many angles
        for jj2=1:repeatNO
            jjTmp=jjTmp+1;
            chosenI=t2(jj2,:,jj1);
            t3{jjTmp}=squeeze(chosenI);

        end
    end

    data3=zeros(o2,n);
    data3Std=data3;        
     for jj=1:o2
        imgRaw=data2(t3{jj},:);
    %         [imgRaw,shiftAll,max2]=imageRegistration(imgRaw,1);
    %         [imgRaw,shiftAll,max2]=imageRegistration_previousRef(imgRaw,1,1);
    %         [imgAll2,shiftAll,max2]=imageRegistration_previousRef(imgAll,upsampling,firstRef)
        img=nanmean(imgRaw,1);
        img2=img;
        data3(jj,:)=img2;
        data3Std(jj,:)=(nanstd(imgRaw,0,1));
        angleLeft2(jj)=angleLeft(t3{jj}(1));



    end      
end

if nargin>5
    if input.SEflag==1
        data3Std=data3Std/sqrt(trialNO);
    end
end

  