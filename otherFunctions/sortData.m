function dff0=sortData(dff0,order,dim)
    
    [m,n]=size(dff0);
%     dff0_sort=zeros(m,n);
    for ii=1:n
        y=dff0(:,ii);
        y=sortData_sub(y,order,dim);
        dff0(:,ii)=y;
    end
    function y=sortData_sub(y,order,dim)
        pp=length(y);
        pAll=order.pAll;
        seq=order.orderAll(:,dim);              
%         dim1=1;
%         for ii=1:dim
%             dim1=dim1*pAll(ii);
%         end
        dim1=prod(pAll(1:dim));
        dim2=pp/dim1;
        y2=reshape(y,dim1,dim2);
        seq2=reshape(seq,dim1,dim2);
        for ii=1:dim2
            seq3=seq2(:,ii);
            [~,I]=sort(seq3);
            y2(:,ii)=y2(I,ii);
        end
        y=y2(:);
        

        