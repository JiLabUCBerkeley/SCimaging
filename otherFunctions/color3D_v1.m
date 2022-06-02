function color3D_v1(filePath,file)
% this function is to convert 3D stack into color coded image. input
% argument filePath is the location where you put file; file is the file
% name of 3dStack; input.evenShift is to correct misalignment between even
% lines and odd lines during bi-directional scanning. 

input.evenShift=0;
if nargin==0
    filePath='C:\serverData\2017.04.14 AM362330 b2gGcamp6f\3D';
    file='*.tif'; 
end
input.weightFlag=0;

fileName=dir(fullfile(filePath,file));
p=length(fileName);
for ii=p:-1:1
    info=imfinfo(fullfile(filePath,fileName(ii).name));
    p1=length(info);
    if p1<20 || p1>500
        fileName(ii)=[];
    end
end
p=length(fileName);
for ii=1:p
    color3D_sub(filePath,fileName(ii).name,input)
end


function color3D_sub(filePath,file,input)
[~,file1,fileType]=fileparts(file);

data1=readSingleTif(fullfile(filePath,file),1,input);
data=double(data1);
[m,n,p]=size(data);
if isfield(input,'weightFlag')
    if input.weightFlag==1
        weight=zeros(p,1);
        x=1:p;x=x(:);
        for ii=1:p
            img=data(:,:,ii);
            weight(ii)=mean(img(:));

            
        end
            pp=polyfit(x,weight,1);pp(1)=pp(1)/2;
            weight=polyval(pp,x);
            weight=weight/max(weight);
            for ii=1:p
                data(:,:,ii)=data(:,:,ii)/weight(ii);
            end
            
    end
end
data=mat2gray(data);
meanImg=nanmean(data,3);
[maxImg,z]=max(data,[],3);
ratio=0.8;
hRaw=linspace(0,ratio,p);

h1=hRaw(z(:));
hImg=reshape(h1,m,n);
sImg=ones(m,n);
meanImg=mat2gray(meanImg);

%% change maxValue to adjust brigtness
maxV=0.55;
meanImg=mat2gray(meanImg,[0,maxV]);
vImg=mat2gray(meanImg);

colorImg=zeros(m,n,3);
colorImg(:,:,1)=hImg;
colorImg(:,:,2)=sImg;
colorImg(:,:,3)=vImg;
meanColor=hsv2rgb(colorImg);
figure(1);clf;imshow(hsv2rgb(colorImg),[])




m=10;n=256;

h=ones(m,1)*linspace(0,ratio,n);
[m,n]=size(h);
s=ones(m,n);
v=ones(m,n)*1;
map2=zeros(m,n,3);
map2(:,:,1)=h;
map2(:,:,2)=s;
map2(:,:,3)=v;
map2=hsv2rgb(map2);
figure(3);clf;imshow(map2,[]);
result='3dColor';
resultPath=fullfile(filePath,result);
if exist(resultPath)==7
else
    mkdir(resultPath)
end
file2=[file1,'_hMax_vMean',fileType];
imwrite(uint8(meanColor*255),fullfile(resultPath,file2),'tif','compression','none');
file2=['map',fileType];
imwrite(uint8(map2*255),fullfile(resultPath,file2),'tif','compression','none');

%% another method: use mean color
% [m,n,p]=size(data);
% data2=zeros(m,n,3,p);
% 
% thresh=0.1;
% for ii=1:p
%     a=ones(m,n)*hRaw(ii);
%     img=data(:,:,ii);
%     a(img<=thresh)=nan;
%     data2(:,:,1,ii)=a;
%     data2(:,:,2,ii)=ones(m,n);
%     data2(:,:,3,ii)=data(:,:,ii);
% end
% data3=nanmean(data2,4);
% data3(:,:,3)=mat2gray(data3(:,:,3));
% 
% data4=hsv2rgb(data3);
% figure(2);clf;imshow(data4,[])
% trash=0;
% file2=['trash_',file1,'_hMean_vMean',fileType];
% imwrite(uint8(data4*255),fullfile(resultPath,file2),'tif','compression','none');
function data=readSingleTif(file,method,input)
if nargin==1
    method=1;
    input.evenShift=0;
elseif nargin==2
    input.evenShift=0;
end

if method==1
    info=imfinfo(file);
p=length(info);
n=info(1).Width;
m=info(1).Height;

data=zeros(m,n,p,'uint16');
for ii=1:p
    img=(imread(file, ii, 'Info', info));
    img(2:2:end,:)=circshift(img(2:2:end,:),[0,input.evenShift]);
    data(:,:,ii)=img;
end
elseif method==2
    InfoImage=imfinfo(file);
    mImage=InfoImage(1).Width;
    nImage=InfoImage(1).Height;
    NumberImages=length(InfoImage);
    data=zeros(nImage,mImage,NumberImages,'uint16');
    FileID = tifflib('open',file,'r');
    rps = tifflib('getField',FileID,Tiff.TagID.RowsPerStrip);

    for i=1:NumberImages
       tifflib('setDirectory',FileID,i);
       % Go through each strip of data.
       rps = min(rps,nImage);
       for r = 1:rps:nImage
          row_inds = r:min(nImage,r+rps-1);
          stripNum = tifflib('computeStrip',FileID,r);
          img = tifflib('readEncodedStrip',FileID,stripNum);
          img(2:2:end,:)=circshift(img(2:2:end,:),[0,input.evenShift]);
          data(row_inds,:,i)=img;
       end
    end
    tifflib('close',FileID);
end
