function Step4_a_getDFF_neuropil_DF(filePath0,file0,options)
%% modifed on 08/12/2016, 
%there is a function in addition to this main code:neuropil.m for neuropil
%subtraction
if nargin==0
%     filePath0='C:\Users\liangy10\Dropbox (HHMI)\superior colliculus manuscript Liang\cortex\KF18_stim25';
%     filePath0='E:\Yajie personal\Lu_testing';
%     filePath0='E:\Yajie personal\Lu_testing\shrink';
%     filePath0='\\dm11\genie\Yajie_GENIE_stuff\project\backupdatafromoldprojects\20180605 417933gad2xai9-scG6s\redrawROI'
    filePath0='C:\Users\yajie.liang\Box\superior colliculus manuscript Liang\data\20190718_singecell_interneurons\m388771_100';
%     filePath0='C:\Users\lur4\Box\superior colliculus manuscript Liang\data\singlecellimaging\m406570_056\testingLu'
    file0='2018*.tif';
    options.WenzhiROI=0;    
    %% choose different methods to calculate the baseline
    bck.method=3; % 1: mode; 2: fixed, e.g., gray period, it uses the first number of options.deleteN to decide the gray period. 3: 10% mimimum;4:to be determined;
    % 5: user input;6:moving average, now it is the average of 30 time
    % points, change at line108
    bck.grayNumber=3;    
    % only for method 5:
%     bck.filePath='D:\Results\20170421 381553 awake\combined_10percentfit70';
%     bck.file='baseline*.mat';
    options.bck=bck;
end
    neuropilH.neuropilMethod=1;% use area between dilationN and dilationN2 as neuropil comtamination. for method 1, it is a square by widthandhight
    %weight decides how much extent is is used to make sure the values
    %remain positive
    options.refreshFlag=1;
%     neuropilH.flag=1;
    neuropilH.dilateN=2;%pixelright neighboring the roi, how many to skip
    neuropilH.width=17;
    neuropilH.height=17;%
    neuropilH.weight=0.7;
    neuropilH.dilateN2=5;
    
    options.neuropilH=neuropilH;
%      chen=neuropilH(avg,bw,xy,dilateN,width,height,dilateN2,neuropilMethod);
%% negative value issue
options.negativeMethod=0; % 1: set negative value as zero when calculating OSI,DSI,GOSI,GDSI; % 2: shift the curve; others: do nothing
%% boundary        
options.LB=[1 -30 5 0.00001 5 -0.2]; %[amp1 theta sigma1 amp2 sigma2 dc];
options.UB=[1.5 360 180 1.5 180 .5]; %[amp1 theta sigma1 amp2 sigma2 dc];
% options.LB=[1 -30 5 0.00001 5 -0.2]; %[amp1 theta sigma1 amp2 sigma2 dc];
% options.UB=[1.5 360 180 1.5 180 .5]; %[amp1 theta sigma1 amp2 sigma2 dc];
%         LB = [1 -30 5 0.00001 5 0]; % origonal value
%         UB = [1.5 360 180 1.5 180 .2];origonal value
fiji_filePath=filePath0; %ROIs drawn in imagej should be saved here with the name"20160707_03*"
filePath{1}=filePath0;
file{1}=file0;
options.filePath2{1}=fiji_filePath;
p=length(filePath);
options.filterMethod.name=0;
for ii=1:p
    if ~isempty(filePath{ii})
        fileName=dir(fullfile(filePath{ii},file{ii}));
        p1=length(fileName);
        for jj=1:p1
            options.jj=jj;
            options.filePath2{jj}=filePath{ii};
            orientationMap_generalized_sub(filePath{ii},fileName(jj),options);
        end
    end
end
function orientationMap_generalized_sub(filePath,fileName,options)
neuropilH=options.neuropilH;
width=neuropilH.width;
height=neuropilH.height;
dilateN=neuropilH.dilateN;

    bck=options.bck;
    if bck.method==5
        bckFileName=dir(fullfile(bck.filePath,bck.file));
        bck.h=matfile(fullfile(bck.filePath,bckFileName(1).name),'writable',true);
        baseline_max=bck.h.baseline_max;
        h.baseline_max=baseline_max;        
    end
disp([filePath,'start----------------------------------------------------------'])
id=1:11;

if options.refreshFlag==1
    method=bck.method;
    data=loadImgSequence(filePath,fileName(1).name);
    try
    time=getTime(fullfile(fileparts(filePath),[fileName(1).name(id),'*.png']));
    catch me
        time=1:size(data,3);
    end
%     data=readSingleTif(fullfile(filePath,fileName(1).name));
    [m,n,p]=size(data);
    [order,qt]=StimulationSequence(filePath,fileName,p);
    filterMethod=options.filterMethod;
    data=differentTypeReadFilter(data,filterMethod);
    order=order{1};
    avg=nanmean(data,3);
    stdImg=nanstd(single(data),0,3);  
    filePathTmp=options.filePath2{options.jj};
    if options.WenzhiROI==1
        [bw,xy]=readWenzhiROI(filePathTmp,fileName(1).name(id),avg);  
    else
        fileName2=dir(fullfile(filePathTmp,[fileName(1).name(id),'*.zip']));

        file2=fileName2(1).name;
        [bw,xy]=readImageJROI_main(filePathTmp,file2,avg);   
        trash=0;
        
    end
%     filePathTmp='C:\Users\lur\Dropbox (HHMI)\manuscript\Bessel beam\figures\fig3_matlab';
%      filePathTmp='C:\Users\lur\Dropbox (HHMI)\manuscript\Bessel beam\figures\fig3_matlab\overlappedROI';

    p1=length(xy);
    Intensity=zeros(p,p1);
    Intensity_subtract=Intensity;
    backIntensity=Intensity;
    dff0=zeros(p,p1);
    f0s=cell(p1,1);
    backF0=cell(p1,1);
    baseline=zeros(p1,1);
    baseline1=zeros(p1,1);
    baseline2=zeros(p,p1);
    tic;
    chen=neuropil(avg,bw,xy,dilateN,width,height,neuropilH.dilateN2,neuropilH.neuropilMethod);
    bwROIAll=chen.bwROIAll;
    neuropilH.baseline=zeros(p1,1);
    backDFF=zeros(p,p1);
    
    I=chen.I;
    for ii=1:p1
          bwROI=bwROIAll(:,:,ii);
        I1=I{ii};
        [Iy,Ix]=ind2sub([m,n],I1);
        Iy2=min(Iy):max(Iy);
        Ix2=min(Ix):max(Ix);
        crop_back=data(Iy2,Ix2,:);
        crop_mask=bwROI(Iy2,Ix2);
        if sum(crop_mask(:))==0
             backIntensity(:,ii)=0;
        else
            backIntensity(:,ii)=squeeze(nanmean(squeeze(nanmean(bsxfun(@times, double(crop_back), double(crop_mask))))))*length(crop_mask(:))/sum(crop_mask(:));
        end
        %% removing F from backIntensity
        
%          if bck.method==1
%             [y2,f0s{ii},baseline(ii),baseline1]=calculatedDff0(Intensity_subtract(:,ii),1);
%         elseif bck.method==2
%             bck.fixPersti=[bck.grayNumber,qt(1)];
%             [y2,f0s{ii},baseline(ii),baseline1]=calculatedDff0(Intensity_subtract(:,ii),bck.method,bck);
%               
%         elseif bck.method==3
%             [y2,f0s{ii},baseline(ii),baseline1]=calculatedDff0(Intensity_subtract(:,ii),4);
%         elseif bck.method==6 %%moving average
%             [y2,f0s{ii},baseline(ii),baseline1]=calculatedDff0(Intensity_subtract(:,ii),6);
%         elseif bck.method==5
%             bck.baseFixed=baseline_max(ii);
%             [y2,f0s{ii},baseline(ii)]=calculatedDff0(Intensity_subtract(:,ii),bck.method,bck);
%         end       
        display(['geting ROI ',num2str(ii)])
        
%         display(['geting ROI ',num2str(ii)])
        
%         speed up the speed
%         y=zeros(p,1);
%         for jj=1:p
%             img=data(:,:,jj);
%             y(jj)=nanmean(img(bw(:,:,ii)));
%         end
        ROIpos = xy{ii};
        crop_frames = data(max(1,floor(min(ROIpos(:,2)))):min(m,ceil(max(ROIpos(:,2)))),max(1,floor(min(ROIpos(:,1)))):min(n,ceil(max(ROIpos(:,1)))),:);
        single_frame_mask = bw(max(1,floor(min(ROIpos(:,2)))):min(m,ceil(max(ROIpos(:,2)))),max(1,floor(min(ROIpos(:,1)))):min(n,ceil(max(ROIpos(:,1)))),ii);
        Intensity(:,ii) = squeeze(nanmean(squeeze(nanmean(bsxfun(@times, double(crop_frames), double(single_frame_mask))))))*length(single_frame_mask(:))/sum(single_frame_mask(:));
        %% neuropil DFF  backDFF
        if bck.method==1
            [y2_neuropil,backF0{ii},neuropilH.baseline(ii)]=calculatedDff0(backIntensity(:,ii),1);
        elseif bck.method==2
            bck.fixPersti=[bck.grayNumber,qt(1)];
            [y2_neuropil,backF0{ii},neuropilH.baseline(ii)]=calculatedDff0(backIntensity(:,ii),bck.method,bck);
              
        elseif bck.method==3
            [y2_neuropil,backF0{ii},neuropilH.baseline(ii)]=calculatedDff0(backIntensity(:,ii),4);
        elseif bck.method==6 %%moving average
            [y2_neuropil,backF0{ii},neuropilH.baseline(ii)]=calculatedDff0(backIntensity(:,ii),6);
        elseif bck.method==5
            bck.baseFixed=baseline_max(ii);
            [y2_neuropil,backF0{ii},neuropilH.baseline(ii)]=calculatedDff0(backIntensity(:,ii),bck.method,bck);
        end
        backDFF(:,ii)=y2_neuropil(:);
        %%
        
        Intensity_subtract(:,ii)=Intensity(:,ii)-backDFF(:,ii)*neuropilH.weight;
        if bck.method==1
            [y2,f0s{ii},baseline(ii),baseline1]=calculatedDff0(Intensity_subtract(:,ii),1);
        elseif bck.method==2
            bck.fixPersti=[bck.grayNumber,qt(1)];
            [y2,f0s{ii},baseline(ii),baseline1]=calculatedDff0(Intensity_subtract(:,ii),bck.method,bck);
              
        elseif bck.method==3
            [y2,f0s{ii},baseline(ii),baseline1]=calculatedDff0(Intensity_subtract(:,ii),4);
        elseif bck.method==6 %%moving average
            [y2,f0s{ii},baseline(ii),baseline1]=calculatedDff0(Intensity_subtract(:,ii),6);
        elseif bck.method==5
            bck.baseFixed=baseline_max(ii);
            [y2,f0s{ii},baseline(ii)]=calculatedDff0(Intensity_subtract(:,ii),bck.method,bck);
        end
        
        dff0(:,ii)=y2(:);
%         if ~isempty(baseline1)%was commented back on feb20, 2019
        baseline2(:,ii)=baseline1(:);
%         end %was commented back on feb20, 2019
    end
    tt=toc;
    weight=neuropilH.weight;
    display(['geting intensities for individual ROIs takes ',num2str(tt,'%.1f'),'seconds'])
    file4=fileName(1).name;
    fileSave=['ROI_fj',strrep(fileName(1).name,'.tif','_Intensity.mat')];
    save(fullfile(filePath,fileSave),'Intensity', 'baseline', 'f0s','baseline2', ...
        'dff0', 'stdImg','avg','bw','xy','p1','qt','order','file4','time','backIntensity','Intensity_subtract','weight');
    
    
else
    fileName2=dir(fullfile(filePath,['ROI_fj',fileName(1).name(id),'*.mat']));
    load(fullfile(filePath,fileName2(1).name));
     weight=neuropilH.weight;
     backDFF=zeros(size(Intensity));
     backF0=cell(size(Intensity,2),1);
  
     
     for ii=1:size(Intensity,2)
         disp(['recalculating dff0 using subtraction weight:',num2str(neuropilH.weight), ' for ROI',num2str(ii)])
          %% neuropil DFF  backDFF
        if bck.method==1
            [y2_neuropil,backF0{ii},neuropilH.baseline(ii)]=calculatedDff0(backIntensity(:,ii),1);
        elseif bck.method==2
            bck.fixPersti=[bck.grayNumber,qt(1)];
            [y2_neuropil,backF0{ii},neuropilH.baseline(ii)]=calculatedDff0(backIntensity(:,ii),bck.method,bck);
              
        elseif bck.method==3
            [y2_neuropil,backF0{ii},neuropilH.baseline(ii)]=calculatedDff0(backIntensity(:,ii),4);
        elseif bck.method==6 %%moving average
            [y2_neuropil,backF0{ii},neuropilH.baseline(ii)]=calculatedDff0(backIntensity(:,ii),6);
        elseif bck.method==5
            bck.baseFixed=baseline_max(ii);
            [y2_neuropil,backF0{ii},neuropilH.baseline(ii)]=calculatedDff0(backIntensity(:,ii),bck.method,bck);
        end
        backDFF(:,ii)=y2_neuropil(:);

         Intensity_subtract(:,ii)=Intensity(:,ii)-backDFF(:,ii)*neuropilH.weight;
          if bck.method==1
            [y2,f0s{ii},baseline(ii),baseline1]=calculatedDff0(Intensity_subtract(:,ii),1);
        elseif bck.method==2
            bck.fixPersti=[bck.grayNumber,qt(1)];
            [y2,f0s{ii},baseline(ii),baseline1]=calculatedDff0(Intensity_subtract(:,ii),bck.method,bck);
        elseif bck.method==3
            [y2,f0s{ii},baseline(ii),baseline1]=calculatedDff0(Intensity_subtract(:,ii),4);
        elseif bck.method==6 %%moving average
            [y2,f0s{ii},baseline(ii),baseline1]=calculatedDff0(Intensity_subtract(:,ii),6);
        elseif bck.method==5
            bck.baseFixed=baseline_max(ii);
            [y2,f0s{ii},baseline(ii)]=calculatedDff0(Intensity_subtract(:,ii),bck.method,bck);
        end
         if mod(ii,5)==1
            t=1:size(Intensity,1);t=t.';
            figure(14);clf;
            subplot(5,1,1);plot(t,backIntensity(:,ii));title(['neuropil, raw, ID= ',num2str(ii)]); hold on; plot(t(backF0{ii}),backIntensity(backF0{ii},ii),'.g');
            subplot(5,1,2);plot(t,y2_neuropil); title('df of neuropil')
            subplot(5,1,3);plot(t,Intensity(:,ii));title('neuron')
            subplot(5,1,4);plot(t,Intensity_subtract(:,ii));title('neuron-df of neuropil')
            hold on;
            plot(t(f0s{ii}),Intensity_subtract(f0s{ii},ii),'.g');
            subplot(5,1,5);plot(t,y2);title('dff of neuron')
            set(gcf,'position',[80         65        1136         813])
            drawnow;
         end
        
        dff0(:,ii)=y2(:);
%         if ~isempty(baseline1)%was commented back on feb20, 2019
        baseline2(:,ii)=baseline1(:);       
     end
         save(fullfile(filePath,fileName2(1).name),'Intensity', 'baseline', 'f0s','baseline2', ...
        'dff0', 'stdImg','avg','bw','xy','p1','qt','order','file4','time','backIntensity','Intensity_subtract','weight','backDFF','neuropilH');
end
