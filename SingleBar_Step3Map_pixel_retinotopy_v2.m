function SingleBar_Step3Map_pixel_retinotopy_v2(filePath,file,options)
    % from raw data, get avg data (averaged across trials); then calculate
    % DF over F. contrast is defined as (biggerHalfMean-smallerHalfMean)/(biggerHalfMean+smallerHalfMean) over all
    % directions. 
if nargin==0

            %%
%     filePath0='C:\serverData\Yajie\20170412 382969sc\s4_singlebar\stack';
    filePath0='G:\superior colliculus manuscript Liang\figures\Figure 1\fig1bc_individualfiles\rawdata';
    file0='20170*.tif';
    options.backgroundIntensity=100;% on final average image    
    options.timeFlag=1;
    options.angleChosen=[2:2:8];     
    options.minContrastFlag=0;options.minContrast=0.6;%     
    options.maxFflag=1; options.maxF=10;% 
    options.dffFlag=1;% use DF over F
    options.deleteN=[0,0];% remove first 2 and last 0
    options.nanROI=[]; % set nan 
    options.refreshFlag=1;  
%     options.minOSI=0.2;
    options.division=3; % divide the stimulation map to 3 by 3;
%      filterMethod_spatial.name=1;filterMethod_spatial.sigmaNumber=3;filterMethod_spatial.sizeNumber=5;
    filterMethod_spatial.name=1;filterMethod_spatial.sigmaNumber=10*2;filterMethod_spatial.sizeNumber=15*2;% 0: do nothing; 1: gaussian; 2: median; 3: wiener    
        filterMethod_temporal.name=1;filterMethod_temporal.sigmaNumber=3;filterMethod_temporal.sizeNumber=5;
% filter; 4: max fitler;5: butterworth; 
    options.filterMethod=filterMethod_spatial;  
    
    options.filterMethod_temporal=filterMethod_temporal;
%     options.WenzhiROI=0;
end
fiji_filePath=filePath0; %ROIs drawn in imagej should be saved here with the name"20160707_03*" 
filePath{1}=filePath0;
file{1}=file0;
options.filePath2{1}=fiji_filePath;
p=length(filePath);
for ii=1:p
    if ~isempty(filePath{ii})
%         matToTifstack_barAdded(filePath{ii},strrep(file{ii},'*_noBar.tif','*.mat'));
        fileName=dir(fullfile(filePath{ii},file{ii}));
        p1=length(fileName);
        for jj=1:p1
            options.jj=jj;
            options.filePath2{jj}=fiji_filePath;
            orientationMap_generalized_sub(filePath{ii},fileName(jj),options);
        end
    end
end

function orientationMap_generalized_sub(filePath,fileName,options)
        options.discreteFlag=1; % always 1
disp([filePath,'start----------------------------------------------------------'])
id=1:11;flagAllH={}; thresh=[0:5:30];
if options.refreshFlag==1
%    method=1;
    data=loadImgSequence(filePath,fileName(1).name);data=single(data);
    avg=mean(data,3);
    filterMethod=options.filterMethod;
    data=differentTypeReadFilter(data,filterMethod); 
    [m,n,p]=size(data);
 if options.dffFlag==1
    
    [m,n,p]=size(data);
    data_avg2=reshape(data,m*n,p);
    data_avg2=data_avg2.';
    data_avg2=calculatedDff0_main(data_avg2,1);
    data=reshape(data_avg2.',m,n,p);    
end   
    flag3=false(m*n,1);
    data_avgTmp=reshape(data,m*n,p);
    
 if options.filterMethod_temporal.name~=0   
    data_avgTmp=data_avgTmp.';
    data_avgTmp=temporalFilter(data_avgTmp,options.filterMethod_temporal);
    data_avgTmp=data_avgTmp.';
    data=reshape(data_avgTmp,m,n,p);
    
 end
     if options.maxFflag==1
         Y_max=max(data_avgTmp,[],2);
        flag3=Y_max<options.maxF;
       
        k=length(thresh);
        
        flagAllH{1}=false(length(flag3),k);
        for ii=1:k
            flagAllH{1}(:,ii)=Y_max<thresh(ii);
        end

%         Y_flag=or(Y_flag,Y_flag_tmp);
    end 
    %% insert nan values
    if ~isempty(options.nanROI)
        data(:,:,options.nanROI)=nan; 
    end;
    [order,point_trial_angle,output]=StimulationSequence(filePath,fileName,size(data,3));    
    trial_angle=point_trial_angle(:,2:3);%  get the trial repetition number and the angle quantity
    framePerSti=point_trial_angle(:,1); %get the the number of frames per stimulus    
    [data_avg,~,angleLeft]=avgTrials(data,order{1}(:,1),framePerSti(1),trial_angle(1,2),trial_angle(1,1));
    file1=fileName(1).name;    
    file2=['bar_',file1];    
%     saveSingleTif(fullfile(filePath,strrep(file2,'bar_','AvgNOBar_')), data_avg)
else
    data_avg=loadImgSequence(filePath,['AvgNOBar_',fileName(1).name]);
    [~,point_trial_angle]=StimulationSequence(filePath,fileName,size(data_avg,3));   
    [order,point_trial_angle]=StimulationSequence(filePath,fileName,point_trial_angle(2)*size(data_avg,3));  
    trial_angle=point_trial_angle(:,2:3);%  get the trial repetition number and the angle quantity
    framePerSti=point_trial_angle(:,1); %get the the number of frames per stimulus        
    [~,~,angleLeft]=avgTrials(order{1}(:,1),order{1}(:,1),framePerSti(1),trial_angle(1,2),trial_angle(1,1));

end


data_avg=single(data_avg);
% avg=mean(data_avg,3);
flag2=avg<=options.backgroundIntensity;
flagAllH{2}=flag2;
 


%% calculate dff0

%% deal with deletion of frames
framePerStiRaw=framePerSti;
data_avg=deleteMN(data_avg,options.deleteN,point_trial_angle(1));
order{1}=deleteMN(order{1},options.deleteN,point_trial_angle(1));
angleLeft=deleteMN(angleLeft,options.deleteN,point_trial_angle(1));
point_trial_angle(:,1)=point_trial_angle(:,1)-sum(options.deleteN);
framePerSti=point_trial_angle(:,1); %get the the number of frames per stimulus 

%%
ver=3;
mapH_p=ceil(framePerSti/ver)*ver;
options.mapH=linspace(0,320,mapH_p)/360;

%% get Map tuned
% data_avgRaw=data_avg;
% [m,n,p]=size(data_avgRaw);
% map=OSI_standard(reshape(data_avgRaw,m*n,p).',point_trial_angle([1,3]),options.minOSI);
% map=reshape(map,m,n);


try
    time=getTime(fullfile(fileparts(filePath),[fileName(1).name(id),'*.png']));
catch me
    time=1:size(data,3);time=time(:);options.timeFlag=0;
end
time=deleteMN(time,options.deleteN,framePerStiRaw);

[time]=avgTrials(time,order{1}(:,1),framePerSti(1),trial_angle(1,2),trial_angle(1,1));

angleNO=point_trial_angle(3);
time=reshape(time,framePerSti,angleNO);
time2=ones(framePerSti,1)*time(1,:);
time=time-time2;
timeAvg=mean(time,2);
if options.timeFlag==1
    staticMovingT=output{1}.staticMovingT;
    timeAvgCrop=timeAvg>staticMovingT(2);
    deleteN=[0,sum(timeAvgCrop)];
    data_avg=deleteMN(data_avg,deleteN,framePerSti);
    framePerSti=framePerSti-sum(deleteN);
end
[m,n,p]=size(data_avg);
data2=reshape(data_avg,m*n,framePerSti,angleNO);
data2b=reshape(data_avg,m*n,framePerSti*angleNO);
[Y,I]=max(data2,[],2);
I=squeeze(I); % 
Y=squeeze(Y); % 
[Y_min,I_min]=min(data2b,[],2);
Y_max=max(data2b,[],2);Y_max=double(Y_max);
if options.minContrastFlag==1
     Y_mean1=zeros(m*n,1);
    Y_mean2=zeros(m*n,1);
    Y_max=0.5*Y_max;
    for ii=1:m*n
        a=data2b(ii,:);
        Y_mean1(ii)=mean(a(a>Y_max(ii)));
        Y_mean2(ii)=mean(a(a<=Y_max(ii)));
    end
    Y_flag=(Y_mean1-Y_mean2)./((Y_mean1+Y_mean2))<options.minContrast;
else
    Y_flag=false(m*n,1);
end
flagAllH{3}=Y_flag;
Y_flag=or(flag3,Y_flag);
% I_min=squeeze(I_min); % 
% Y_min=squeeze(Y_min); % 
% Y_min=double(Y_min);
% Y=double(Y);
% %% 
% Y_tmp=Y_min+Y_max;
% Y_tmp(Y_tmp==0)=1;
% Y_OSI=Y_max-Y_min;
% Y_OSI=Y_OSI./Y_tmp;
% minContrast=options.minContrast;
% Y_flag=Y_OSI<minContrast;
% I(~Y_flag)=nan;
% s=4;


for s=3:5
    if mod(angleNO,s)==0
        I2b=reshape(I,m,s*n,angleNO/s);
        I2=zeros(m*angleNO/s,s*n);
        for ii=1:angleNO/s
            iia=(1:m)+(ii-1)*m;
            I2(iia,:)=I2b(:,:,ii);
        end
        break;
    end    
end
% figure(1);clf;imshow(I2,[])
if filterMethod.name==0
    result='rtnpyMap';
else
    result=['rtpyMap_flted',num2str(filterMethod.name)];
end
result=[result,'dim',num2str(options.division),'Angle',num2str(length(options.angleChosen))];
% options.chosen
resultPath=fullfile(filePath,result);

if exist(resultPath)
else
    mkdir(resultPath)
end
I_save=uint8(reshape(I,m,n,angleNO));
file=fileName(1).name;
saveSingleTif(fullfile(resultPath,strrep(file,'.tif','peakTime.tif')),I_save);

[rgbImg,map2,~,discreteMap]=mapHSV(zeros(size(I2))+1,I2,options);
k1=6;k2=4*4;
kk=1:k1*k2;kk=reshape(kk,k1,k2);
kk2=kk(1:k1-2,:);
figure(2);clf;imshow(uint8(rgbImg*255))
imwrite(uint8(rgbImg*255),fullfile(resultPath,strrep(file,'.tif','peakTime_color.tif')),'tif','compression','none');
figure(3);clf;imshow(uint8(discreteMap*255))
imwrite(uint8(discreteMap*255),fullfile(resultPath,strrep(file,'.tif','peakTime_colorMap.tif')),'tif','compression','none');
yx=zeros(m*n,2);
A2=(1:angleNO)-1;A2=A2*360/angleNO;
A2=A2/360*2*pi;
A=[sin(A2(:)),cos(A2(:))];
r=(I-0.5*framePerSti)/framePerSti;
angleChosen=options.angleChosen;
A=A(angleChosen,:);
r=r(:,angleChosen);

for ii=1:m*n
    if mod(ii,m*5)==0
        disp(['analyzing vertical line',num2str(ii/m,'%03d')])
    end
    B2=r(ii,:);
    B=B2(:);
    X=A\B;
    yx(ii,:)=X.';
    
end
yx(yx>0.5)=0.5;
yx(yx<-0.5)=-0.5;
yx=yx+0.5;

k=options.division;
yx=yx*k;
yx=yx+0.49;
yx=round(yx);
y=yx(:,1);y(y==0)=1;
x=yx(:,2);x(x==0)=1;
ind=sub2ind([k,k],y,x);
indRaw=ind;
% Y_flag2=Y_flag;
% Y_flag3=zeros(size(Y_flag2,1));
% Y_flag2=(sum(Y_flag,2)==8);
% Y_flag4=isnan(Y_flag3);
% Y_flag2=reshape(Y_flag3,m,n,angleNO);
ind(Y_flag)=nan;
ind(flag2(:))=nan;
ind2=reshape(ind,m,n);
options.mapH=([43,233,177,311,58,118,1,280,80]-1)/360;
if k>=4
    options.mapH=linspace(0,320,k*k)/360;
end
[rgbImg,map2,discreteMap]=mapHSV(ind2*0+1,ind2,options);
figure(4);clf;imshow(uint8(rgbImg*255))
imwrite(uint8(rgbImg*255),fullfile(resultPath,['retinpy_',file]),'tif','compression','none');
figure(5);clf;imshow(uint8(discreteMap*255))
imwrite(uint8(discreteMap*255),fullfile(resultPath,['retinpy_map',file]),'tif','compression','none');
trash=0;
k=size(flagAllH{1},2);
p=length(flagAllH);
for ii=1:k
    flag=flagAllH{1}(:,ii);
    for jj=2:p
        flag=or(flag,flagAllH{jj}(:));

    end
    
        ind=indRaw;
        ind(flag)=nan;
        ind2=reshape(ind,m,n);    
        [rgbImg,map2,discreteMap]=mapHSV(ind2*0+1,ind2,options);
        imwrite(uint8(rgbImg*255),fullfile(resultPath,['retinpy_maxF',num2str(thresh(ii),'%02d'),'_',file]),'tif','compression','none');
end



