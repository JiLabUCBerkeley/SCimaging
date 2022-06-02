function data=differentTypeReadFilter(data,filterMethod)
[~,~,p]=size(data);
for ii=1:p
    data(:,:,ii)=differentTypeReadFilterIndividual(data(:,:,ii),filterMethod);
end
function img=differentTypeReadFilterIndividual(img,filterMethod)
% filterMethod.name, 0: do nothing; 1: gaussian; 2: median; 3: wiener
% filter; 4: max fitler;5: butterworth; 
% filterMethod.name=1;filterMethod.sigmaNumber=2;filterMethod.sizeNumber=3;
name=filterMethod.name;


if name<=4 && name>0
    sizeNumber=filterMethod.sizeNumber;
end
if name==5
   order=filterMethod.order;
   fcut=filterMethod.fcut;    
end
%% filterImg
if name==0
    %imgAvg=imgAvg;
elseif name==1
    w=fspecial('gaussian',[filterMethod.sizeNumber,filterMethod.sizeNumber],filterMethod.sigmaNumber);
    img=imfilter(img(:,:,1),w,'replicate');
elseif name==2
    img=medfilt2(img(:,:,1),[sizeNumber,sizeNumber]);
elseif name==3
        img=wiener2(img(:,:,1),[sizeNumber,sizeNumber]);
elseif name==4
    imgAvg2=maxFilter2(img(:,:,1),sizeNumber);
    img=single(imgAvg2);
%     imgAvg2=edge(img(:,:,1),'canny',[.05,.4],sigmaNo);
%     img=single(imgAvg2);
elseif name==5

   img=butter_fcut(img,fcut,1,order,'abs');
end

function [zoFilter,varargout]=butter_fcut(z0,fcut,method,butterN,rc)
% method=1 butterworth;
[m,n]=size(z0);
mm=1:m;
nn=1:n;
m0=m/2;
if mod(m0,2)==0
    m0=m0+0.5;
else
    m0=floor(m0)+1;
end

n0=n/2;
if mod(n0,2)==0
    n0=n0+0.5;
else
    n0=floor(n0)+1;
end
[xx,yy]=meshgrid(nn,mm);
yy=yy-m0;
xx=xx-n0;
dist=sqrt(yy.^2+xx.^2);
if method==1

    map=1 ./ (1.0 + (dist ./ (fcut*m)).^(2*butterN));
    zoFilter=ifft2(ifftshift(map.*fftshift(fft2(z0))));
elseif method==2
    map=dist<fcut*m;
    zoFilter=ifft2(ifftshift(map.*fftshift(fft2(z0))));    
end
if strcmp(rc,'abs')
    zoFilter=abs(zoFilter);
end
if nargout>=2
    varargout{1}=map;
end

function  img2=maxFilter2(img1,sizNo)
[mm,nn]=size(img1);
sizNo=floor(sizNo/2);
edgeY=sizNo;
edgeX=sizNo;
%% make the images bigger; replicate padding
% img1_big=zeros(mm+2*edgeY,nn+2*edgeX);
% img1_big(edgeY+1:edgeY+mm,edgeX+1:edgeX+nn)=img1;
img1_big=zeros(mm+2*edgeY,nn+2*edgeX);
img1_big(edgeY+1:edgeY+mm,edgeX+1:edgeX+nn)=img1;
img1_big(1:edgeY,:)=ones(edgeY,1)*img1_big(1+edgeY,:);
img1_big(edgeY+mm+1:end,:)=ones(edgeY,1)*img1_big(mm,:);
img1_big(:,1:edgeX)=img1_big(:,edgeX+1)*ones(1,edgeX);
img1_big(:,nn+1+edgeX:end)=img1_big(:,nn)*ones(1,edgeX);
sizeY1=edgeY;
sizeX1=edgeX;
%% initialize matrix
img2=zeros(mm,nn);
for iiTmp=1:mm
%          waitbar(iiTmp/mm,h_wait,[num2str(100*iiTmp/mm,'%04.1f'),'%completed']);
    ii=edgeY+iiTmp;
    for jjTmp=1:nn
        jj=edgeX+jjTmp;
        ROI=img1_big(ii-sizeY1:ii+sizeY1,jj-sizeX1:jj+sizeX1); %#ok<*PFBNS>
        img2(iiTmp,jjTmp)=max(ROI(:));
%         disp(['y',num2str(iiTmp,'%03d'),'jj',num2str(jjTmp,'%03d')])
    end
end