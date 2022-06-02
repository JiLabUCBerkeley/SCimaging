function data3D=groupZprojected(data,avgNumber)
[m,n,p]=size(data);
k1=floor(p/avgNumber);
k2=mod(p,avgNumber);
data3D=zeros(m,n,k1);
if k2==0
    i1=cell(k1,1);
    for ii=1:k1
        i1{ii}=[(ii-1)*avgNumber+1:(ii-1)*avgNumber+avgNumber];
    end
else
    i1=cell(k1+1,1);
    for ii=1:k1
        i1{ii}=[(ii-1)*avgNumber+1:(ii-1)*avgNumber+avgNumber];
    end   
    i1{k1+1}=[(k1+1-1)*avgNumber+1:p];
end
k1=length(i1);
for ii=1:k1
    x=data(:,:,i1{ii});
    data3D(:,:,ii)=mean(x,3);

end