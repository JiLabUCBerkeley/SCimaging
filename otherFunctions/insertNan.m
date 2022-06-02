function data2=insertNan(data,repeatNO,ratio)
[m,n]=size(data);
k=round(repeatNO*ratio);
p=m*(k+repeatNO)/repeatNO;
data2=nan(p,n);
allSti=floor(m/repeatNO+eps);
nanFlag=[ones(k,1);zeros(repeatNO,1)]*ones(1,allSti);
nanFlag2=nanFlag(:);

for ii=1:n
    data2(nanFlag2==0,ii)=data(:,ii);
end
