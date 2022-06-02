function test1=movingdeleteMN(test,deleteN,repeatNO)

[m,n,p]=size(test);
if p==1
    allSti=floor(m/repeatNO);
    a=zeros(repeatNO,1);
    a((deleteN(1)+1):end)=1; %%delete visual evoked frame
    deleteFlag=a*ones(1,allSti);
    deleteFlag2=deleteFlag(:);
    test(deleteFlag2==1,:)=[];
       
else
    m=p;
    allSti=floor(m/repeatNO);
    a=zeros(repeatNO,1);
    a(1:deleteN(1))=1;
    if deleteN(2)>0
        a(end-deleteN(2)+1:end)=1;
    end

    deleteFlag=a*ones(1,allSti);
    deleteFlag2=deleteFlag(:);
    test(:,:,deleteFlag2==1)=[];    
end
test1=test;