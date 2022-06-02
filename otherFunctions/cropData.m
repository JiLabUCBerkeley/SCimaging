function dff0_crop=cropData(dff0,deleteN,pAll)
[m,n]=size(dff0);
for ii=1:n
    y=dff0(:,ii);
    y2=deleteMN(y,deleteN,pAll);
    if ii==1
        dff0_crop=zeros(length(y2),n);
    end
    dff0_crop(:,ii)=y2(:);
end