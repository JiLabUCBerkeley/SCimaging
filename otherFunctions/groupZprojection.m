function data3D=groupZprojection(imgRaw,avgNumber, options)
%% to be finished
if nargin==2
    options.registerFlag=0;
end
p=size(imgRaw,3);
if avgNumber==0
   
    data3D=imgRaw;
    return;
    
end


k1=floor(p/avgNumber);
k2=mod(p,avgNumber);
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
[mm,nn]=size(imgRaw(:,:,1));
data3D=zeros(mm,nn,k1);
for ii=1:k1
    if length(i1{ii})==1 || input.registerFlag==0
        x=imgRaw(:,:,i1{ii});
    else        
        [x]=imgRegistration(x);
    end 
    imgMean=mean(x,3);
    data3D(:,:,ii)=imgMean;
end