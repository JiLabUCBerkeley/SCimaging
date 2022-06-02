function [data3,data3Std,angleLeft2]=avgRepeatNO(data2,angleLeft,input)
[m,n,p]=size(data2);
if p>1 % dealing wiht images
    c=unique(angleLeft);
    p2=length(c);

    data3=zeros(m,n,p2,'uint16');
    data3Std=zeros(m,n,p2);
    for ii=1:p2
        cI=(angleLeft==c(ii));
        img=nanmean(double(data2(:,:,cI)),3);
        img2=nanstd(double(data2(:,:,cI)),0,3);
        data3(:,:,ii)=uint16(img);
        data3Std(:,:,ii)=img2;
        if nargin>2
            if input.SEflag==1
                data3Std(:,:,ii)=data3Std(:,:,ii)/sqrt(sum(cI));% I used length(cI) before, which was severely wrong. 
            end
        end        
    end 
    angleLeft2=c;    
else % dealing with curves
    c=unique(angleLeft);
    p2=length(c);
    data3=zeros(p2,n);
    data3Std=zeros(p2,n);
    for ii=1:p2
        cI=(angleLeft==c(ii));
        img=nanmean(double(data2(cI,:)),1);
        img2=nanstd(double(data2(cI,:)),0,1);
        data3(ii,:)=img;
        data3Std(ii,:)=img2;
        if nargin>2
            if input.SEflag==1
                data3Std(ii,:)=data3Std(ii,:)/sqrt(sum(cI));
            end
        end         
    end
    angleLeft2=c; 
    
end

