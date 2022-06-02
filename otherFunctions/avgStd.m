function [avg,SE]=avgStd(dff0,pAll,dim)
    % calculate average and std along multiple dimensions; pAll indicates
    % the length for each dimensions; dim is a vector, containing the
    % dimensions for average and SE calculation. 
    [~,n]=size(dff0);
            pAll2=1:length(pAll);
            pAll2(dim)=[];    
    for ii=1:n
        y=(dff0(:,ii));
        y2=reshape(y,pAll);
        if length(dim)==1
            avgY=mean(y2,dim);
            SEY=std(y2,0,dim);

        else

            y3=permute(y2,[dim(:).',pAll2(:).']);
            y4=reshape(y3,prod(pAll(dim)),prod(pAll(pAll2)));
            avgY=mean(y4,1);
            SEY=std(y4,0,1);

        end
%         SEY2=SEY(:);
        avgY2=avgY(:);
        SEY2=SEY(:);
        SEY2=SEY2/sqrt(prod(pAll(dim)));
        if ii==1
            avg=zeros(length(avgY2),n);
            SE=avg;
            
        end
        avg(:,ii)=avgY2;
        SE(:,ii)=SEY2;
        
    end
