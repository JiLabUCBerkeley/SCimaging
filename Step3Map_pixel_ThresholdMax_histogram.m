function Step3Map_pixel_ThresholdMax_histogram(filePath0,file0,options)
% generate pixel maps for orientations; it need two inputs:the matlab file containing visual stimulation setting and the raw tif
if nargin==0
%         filePath0='\\dm11\genie\Yajie_GENIE_stuff\project\backupdatafromoldprojects\20180530 401905ntsrGN209\02\stack'; 
%     file0='2018*.tif';
    filePath0='G:\superior colliculus manuscript Liang\figures\Figure 1\fig1d_rawfiles\rawdata_generatingOSmap';
    file0='2018*.tif';
    filterMethod.name=1;% 0: do nothing; 1: gaussian; 2: median; 3: wiener
    % filter; 4: max fitler;5: butterworth; 
    filterMethod.name=1;filterMethod.sigmaNumber=1;filterMethod.sizeNumber=3;
    options.filterMethod=filterMethod;
    options.deleteN=[3,0];%[18,0];% remove first 2 and last 0
    options.nanROI=nan; % set nan 
     bck.method=2; % 1: mode; 2: fixed, e.g., gray period, it uses the first number of options.deleteN to decide the gray period. 3: 10% mimimum;4:to be determined;
     bck.grayNumber=options.deleteN(1);
     options.bck=bck;
%% negative value issuep
options.negativeMethod=0; % 1: set negative value as zero when calculating OSI,DSI,GOSI,GDSI; % 2: shift the curve; others: do nothing
%% boundary
options.LB=[1 -30 5 0.00001 5 -0.2]; %[amp1 theta sigma1 amp2 sigma2 dc];
options.UB=[1.5 360 180 1.5 180 .5]; %[amp1 theta sigma1 amp2 sigma2 dc];
%         LB = [1 -30 5 0.00001 5 0]; % origonal value
%         UB = [1.5 360 180 1.5 180 .2];origonal value    

end
criteria.maxYThreshFlag=1; 
criteria.maxYThresh=10; % threshold to determine weather the cell is responsive; default value: 10 (10%)        
options.criteria=criteria;

filePath{1}=filePath0;
file{1}=file0;

p=length(filePath);
for ii=1:p
    if ~isempty(filePath{ii})
        fileName=dir(fullfile(filePath{ii},file{ii}));
        p1=length(fileName);
        for jj=1:p1
            orientationMap_generalized_sub(filePath{ii},fileName(jj),options);
        end
    end
end

function orientationMap_generalized_sub(filePath,fileName,input)
criteria=input.criteria;
data=loadImgSequence(filePath,fileName(1).name);
% data=readSingleTif(fullfile(filePath,fileName(1).name));
[m,n,p]=size(data);
[order,point_trial_angle]=StimulationSequence(filePath,fileName,p);
filterMethod=input.filterMethod;
data=differentTypeReadFilter(data,filterMethod);
order=order{1};
input.kk=1;%% spines
avg=nanmean(data,3);
stdImg=nanstd(single(data),0,3);
%  [a, MSGID] = lastwarn();
warning('off','stats:statrobustfit:IterationLimit');
bw1=true(m,n);
% method 1, mod; most frequent mod
%method 2 fixed: fixPersti has two inputs; fixPersti(1) is the static frame number per stimulus;
% fixPersti(2) is the moving frame number per stimulus;
% method 3: fix is a vector, defing the frame numbers as baseline; example:
% fix=[1,10,11,12], then average the frame # 1, 10, 11, 12 as baseline
% method 4:  use the minimum 20% values as baseline;
method=1;%mod
deleteN=input.deleteN;
nanROI=input.nanROI;
angleNO=point_trial_angle(3);
trialNO=point_trial_angle(2);
framePerSti2=length(order)/angleNO/trialNO;
order=deleteMN(order,deleteN,framePerSti2);
[Y,I]=sort(order(:,1));
order(:,2:4)=[(1:length(I)).',Y(:)*1,I(:)]; 
framePerSti=length(order)/angleNO/trialNO;
point_trial_angle=[framePerSti,angleNO,trialNO];
input.angleNO=angleNO;
input.trialNO=trialNO;
input.framePerSti=framePerSti;
input.qt=point_trial_angle;
pvalueMin=0.05/1;
map=nan(m,n);
pvalue=ones(m,n);
tic;
poolobj = gcp('nocreate');
if isempty(poolobj)
    c=parcluster;
    parpool('local',c.NumWorkers);
else
%     poolsize = poolobj.NumWorkers
end
% for ii=1:m
bck=input.bck;
if bck.method==5
    bckFileName=dir(fullfile(bck.filePath,bck.file));
    bck.h=matfile(fullfile(bck.filePath,bckFileName(1).name),'writable',true);
    baseline_max=bck.h.baseline_max;
    h.baseline_max=baseline_max;
else
    baseline_max=[];

end
% for ii=1:m    
parfor ii=1:m
    bck2=bck;
    disp(['analyzing line',num2str(ii,'%03d')])
    for jj=1:n
        if bw1(ii,jj)
            % obtain temporal curve for individual pixels;
            y=squeeze(data(ii,jj,:));% reduce the redundant dimension
            y2=double(y(:));
%             y2=calculatedDff0(y2,method);% change to dff0
        if bck2.method==1
            y2=calculatedDff0(y2,method);% change to dff0
        elseif bck2.method==2
            bck2.fixPersti=[bck2.grayNumber,framePerSti2];
            y2=calculatedDff0(y2,bck2.method,bck2);% change to dff0
        elseif bck2.method==3
            y2=calculatedDff0(y2,4);%
        elseif bck2.method==5
%             bck2.baseFixed=baseline_max(ii);
%             y2=calculatedDff0(y2,bck2.method,bck2);
        end
%             brob=robustfit(imgM,y2);

%             y2=y2-(brob(1)+brob(2)*imgM);

            if length(nanROI)>1
                y2(nanROI,:)=nan;
            else
            end
            
            y3=deleteMN(y2,deleteN,framePerSti2);
%             y3=;
% get pvalue for anova test
            pvalue(ii,jj)=pvalueGet(y3(order(:,end)),input); 
            % do double Gaussian fitting if pvalue is small than pvalueMin
            if pvalue(ii,jj)<=pvalueMin
                output=GaussianFit(y3,order,input);
                
                 if criteria.maxYThreshFlag==1
                    maxY=output.maxY;
                    if  max(maxY)>=criteria.maxYThresh
                        map(ii,jj)=output.theta;
                    end
                 else
                     map(ii,jj)=output.theta;
                end               
            end
        end
        


        
    end
end
%% saveMap
    fileSave=['histogram',strrep(fileName(1).name,'.tif','.mat')];
    save(fullfile(filePath,fileSave),'map', 'input', 'pvalueMin','order');
%% processing avg
avg=mat2gray(avg);
avgTmp=sort(avg(:));
ratioAvg=0.99;
avgUpperLimit=avgTmp(round(ratioAvg*m*n));
avg=mat2gray(avg,[0,avgUpperLimit]);

[rgbImg,map2]=mapHSV(avg,map);
stdImg2=mapHSV(stdImg,map);
if filterMethod.name==0
    result='pixelMap';
else
    result=['pixelMap_flted',num2str(filterMethod.name)];
end
if input.deleteN(1)>0
    result=[result,'_trimp',num2str(input.deleteN(1))];
end
if criteria.maxYThreshFlag==1
    result=[result,'_max',num2str(criteria.maxYThresh)];
end
resultPath=fullfile(filePath,result);
if exist(resultPath)
else
    mkdir(resultPath)
end
id=1:11;
file=[fileName(1).name(id),'_colormap.tif'];
cmbFile=['cmb_',fileName(1).name(id),'.tif'];
imwrite(uint8(rgbImg*255),fullfile(resultPath,['cmb_',file]),'tif','compression','none');
file=fullfile(resultPath,file);
imwrite(uint8(stdImg2*255),strrep(file,'.tif','_std.tif'),'tif','compression','none');

file1=strrep(file,'map.tif','pValue_color.mat');
save(file1,'map','pvalue','rgbImg','stdImg2');
%%
img=avg*0+255;
img=uint8(img);

[rgbImg,map2]=mapHSV(img,map);
[mm,nn,p]=size(rgbImg);
imgAll=zeros(mm,nn*3+4,p);
% imgAll(:,:,1)=0.5;
s=0;
s=s+1;imgAll(:,(s-1)*nn+(1:nn)+s-1,:)=rgbImg;
figure(2);clf;subplot(1,3,1);
imshow(uint8(rgbImg*256));title('360')
figure(3);subplot(1,3,1);imshow(map2)
imwrite(uint8(rgbImg*255),strrep(file,'.tif','_360Fullmap.tif'),'tif','compression','none');
%%
[rgbImg,map2]=mapHSV_180Degrees(img,map);
figure(2);
subplot(1,3,2);
imshow(uint8(rgbImg*256));title('180');
[m,n,~]=size(map2);
figure(3);subplot(1,3,2);imshow(map2(:,1:round(n/2),:))

s=s+1;imgAll(:,(s-1)*nn+(1:nn)+s,:)=rgbImg;

imwrite(uint8(rgbImg*255),strrep(file,'.tif','_180Halfmap.tif'),'tif','compression','none');
[rgbImg,map2]=mapHSV_180Degrees(avg,map);
imwrite(uint8(rgbImg*255),strrep(fullfile(resultPath,cmbFile),'.tif','_180Halfmap.tif'),'tif','compression','none');
%%
input.fullColor=1;
input.lhFlag=0;
[rgbImg,map2]=mapHSV_180Degrees(img,map,input);
imwrite(uint8(rgbImg*255),strrep(file,'.tif','_180Fullmap.tif'),'tif','compression','none');
% combined raw figure and pixel map with 180 degree full colormap
[rgbImg2,map2]=mapHSV_180Degrees(avg,map,input);
imwrite(uint8(rgbImg2*255),strrep(fullfile(resultPath,cmbFile),'.tif','_180Fullmap.tif'),'tif','compression','none');

s=s+1;imgAll(:,(s-1)*nn+(1:nn)+s,:)=rgbImg;
figure(2);
subplot(1,3,3);
imshow(uint8(rgbImg*256));title('180');
[m,n]=size(map2);
figure(3);subplot(1,3,3);imshow(map2);
file5=strrep(file,'.tif','_blank.eps');
figure(2);
print('-r300',strrep(file5,'.eps','.png'),'-dpng')
set(gcf,'PaperPositionMode','auto')
print('-r300',file5,'-depsc','-tiff')
figure(3);
file5=strrep(file,'.tif','_clcbar.eps');
print('-r300',strrep(file5,'.eps','.png'),'-dpng');
set(gcf,'PaperPositionMode','auto')
print('-r300',file5,'-depsc','-tiff')
imwrite(uint8(imgAll*255),strrep(file,'.tif','_blk_img.tif'),'tif','compression','none');
function shadedErrorStack(data)
[m,n]=size(data);
data2=data;
roiN=(n-1)/2;% first column is x;
x=data(:,1);
s=0;
step=zeros(roiN,1);
ratio=ones(roiN,1);
for ii=2:roiN
    maxV=max(max(data(2:end)));
    step(ii)=maxV*1;
end
maxV2=zeros(roiN,1);
input.linepros=cell(roiN,1);
for ii=1:roiN
    dataIn=[x,data(:,2*ii:2*ii+1)*ratio(ii)];
    maxV2(ii)=max(abs(dataIn(:,2)));
    input.linepros{ii}={'-k','lineWidth',0.7};
end
for ii=1:2:roiN
    input.linepros{ii}={'-r','lineWidth',0.7};
end
maxV2(2:2:end)=maxV2(1:2:end);
maxV2=maxV2/max(maxV2);
ratio=1./maxV2;
ratio(ratio<=2)=1;
ratio(ratio>2&ratio<=4)=2;
ratio(ratio>4)=3;
% ratio=ones(roiN,1);
fig=11;
figure(fig);clf;
input.fig=fig; hold on;
 xx=1;
 scaleB=50;
for ii=1:roiN
    dataIn=[x,data(:,2*ii:2*ii+1)*ratio(ii)];
    dataIn(:,2)=dataIn(:,2)-sum(step(1:ii));
    shadedErrorStack_sub(dataIn,input.linepros{ii});
    
    y0=nanmean(dataIn(:,2));
    scaleB2=scaleB*ratio(ii);
    y1=y0-scaleB2/2;
    y2=y1+scaleB2-1;
    plot([1,1]*xx,[y1,y2],input.linepros{ii}{:})
end
ylim([-5200,200])
set(gcf,'PaperPositionMode','auto')
set(gcf,'position',[553   167   612   946]);
 set(gca,'TickDir','out');
 axis off;

%  for ii=1:roiN
%      plot([1,1]*x,[50,y0+scaleB*ratio(ii)]-sum(step(1:ii)),'k-')
%  end
%  plot([1,1],[100,200],'k-')
%    
function shadedErrorStack_sub(dataIn,linepros)
x=dataIn(:,1);
p=length(x);
y=dataIn(:,2);
ySE=dataIn(:,3);
b=~isnan(ySE);
c=b(2:end)-b(1:end-1);
x1=find(c==1);x1=x1+1;
x2=find(c==-1);
x1=x1(:);
x2=x2(:);
if b(1)
    x1=[1;x1];
end
if b(end)
    x2=[x2;p];
end
xx=[x1,x2];
p1=size(xx,1);
% linepros={'-k','lineWidth',2};
for ii=1:p1
    xIn=xx(ii,1):xx(ii,2);
    shadedErrorBar(x(xIn),y(xIn),ySE(xIn),linepros);
end
trash=0;

function dataOut=dataExtract(order,aa,dataFit,input)
nanRatio=1;
deleteN=[2,0];
angleNO=12;
trialNO=10;
repeatNO2=length(order)/angleNO/trialNO;
order=deleteMN(order,deleteN,repeatNO2);
[Y,I]=sort(order(:,1));
order(:,2:4)=[(1:length(I)).',Y(:)*1,I(:)]; 
repeatNO=length(order)/angleNO/trialNO;
nanROI=[900:906];
dataOut=cell(length(dataFit.branch),4);
for s=1:length(dataFit.branch)
    test=dataFit.branch{s}.dff0(:,aa{s}(:,input.kk)); 
    test(nanROI,:)=nan;

    % test=dataFit.branch{2}.dff0(:,b1(:,1));
    % test=dataFit.branch{3}.dff0(:,c1(:,1));
    % test=dataFit.branch{4}.dff0(:,d1(:,1));


    
    %% delete first first 2 time points
    test=deleteMN(test,deleteN,repeatNO2);
    
%     
     input.SEflag=1;
    [data3,data3SE,angleLeft2]=avgTrials(test,order(:,1),repeatNO,angleNO,trialNO,input);
    data3=insertNan(data3,repeatNO,nanRatio);
    data3SE=insertNan(data3SE,repeatNO,nanRatio);
    dataOut{s,1}=[(1:length(data3)).',data3,data3SE];
    dataOut{s,1}(:,2:2:end)=data3;
    dataOut{s,1}(:,3:2:end)=data3SE;
    x=1:size(data3,1);
    [data3b,data3SEb,angleLeft2b]=avgRepeatNO(test(order(:,end),:),order(:,3),input);
    
    x2=1:size(data3b,1);x2=x2-1;x2=x2*30;
    dataOut{s,2}=[x2(:),data3b,data3SEb];
    dataOut{s,2}(:,2:2:end)=data3b;
    dataOut{s,2}(:,3:2:end)=data3SEb;
    ii=1;
    xx=1:size(data3b,1);xx=xx-1;xx=xx*360/size(data3b,1);
    for ii=1:size(test,2)
%         figure(1);clf; subplot(2,1,1);
% %         plot(x,data3(:,ii))
%         errorbar(x,data3(:,ii),data3SE(:,ii));
%         subplot(2,1,2);
%         errorbar(x2,data3b(:,ii),data3SEb(:,ii),'ob');
% 
%         hold on;
        normFlag=1;
        output=curveGaussianFitting(xx,data3b(:,ii),data3SEb(:,ii),normFlag);
        if ii==1
            dataOut{s,3}=zeros(length(output.x2),1+size(test,2));
            dataOut{s,3}(:,1)=output.x2(:);
            dataOut{s,4}=zeros(size(test,2),1);
        end
        dataOut{s,4}(ii)=output.theta;
        dataOut{s,3}(:,ii+1)=output.z2(:);
%         plot(output.x2,output.z2);
%         ylim0=get(gca,'yLim');
%         plot([1 1]*output.theta,ylim0,'r-')
%         title(['branch=',num2str(s),';id=',num2str(ii,'%02d'),';spine=',num2str(aa{s}(ii,1)),'angle=',num2str(round(output.theta))]);
    end
end



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


function [data3,data3SE]=trialSortAvg(data,s,order)
data2=data;
[m,n]=size(data);

for ii=1:n
    data2(:,ii)=data(order(:,end),ii);
end
p2=floor(m/s);
s2=ones(m,1)*p2;

s3=(1:s).'*ones(1,p2);
s3b=s3(:);
p3=length(s3b);
s2(1:p3)=s3b(:);
input.SEflag=1;
[data3,data3SE,angleLeft2]=avgRepeatNO(data2,s2,input);


function b1=findDentrite(a,a1)
[m,~]=size(a);
a1=a1(:);
[m2,~]=size(a1);
b1=zeros(m2,2);
b1(:,1)=a1;
for ii=1:m2
    s=find(a(:,1)==a1(ii));
    b1(ii,2)=a(s,2);
end


