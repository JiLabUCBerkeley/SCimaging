function Step4_a_getDFF(filePath0,file0,fiji_filePath,options)
%% modifed on 08/12/2016, 
% 1) add an option to have Wenzhi's ROIs
% 2) create tuning map by putting color of tuning orientations for individual ROIs
% 3) add ROI number on top of tuning map
% 4) add tuning angles on top of tuning map
% modifed on 08/15/2015,
% correct for the batch processing error;
if nargin==0
    filePath0='C:\Users\yajie.liang\Box\superior colliculus manuscript Liang\figures\Figure 1\rawfiles\tuningmapandcurve';
    file0='2018*.tif';
    options.WenzhiROI=0;    
    %% choose different methods to calculate the baseline
    bck.method=2; % 1: mode; 2: fixed, e.g., gray period, it uses the first number of options.deleteN to decide the gray period. 3: 10% mimimum;4:to be determined;
    % 5: user input
    bck.grayNumber=3;    
    % only for method 5:
    bck.filePath='D:\Results\20170421 381553 awake\combined_10percentfit70';
    bck.file='baseline*.mat';
    options.bck=bck;
end
%% negative value issue
options.negativeMethod=2; % 1: set negative value as zero when calculating OSI,DSI,GOSI,GDSI; % 2: shift the curve; others: do nothing
%% boundary        
options.LB=[1 -30 5 0.00001 5 -0.2]; %[amp1 theta sigma1 amp2 sigma2 dc];
options.UB=[1.5 360 180 1.5 180 .5]; %[amp1 theta sigma1 amp2 sigma2 dc];
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
    bck=options.bck;
    if bck.method==5
        bckFileName=dir(fullfile(bck.filePath,bck.file));
        bck.h=matfile(fullfile(bck.filePath,bckFileName(1).name),'writable',true);
        baseline_max=bck.h.baseline_max;
        h.baseline_max=baseline_max;        
    end
disp([filePath,'start----------------------------------------------------------'])
id=1:11;
options.refreshFlag=1;
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
    dff0=zeros(p,p1);
    f0s=cell(p1,1);
    baseline=zeros(1,p1);
    tic;
    for ii=1:p1
        display(['geting ROI ',num2str(ii)])
        
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
        if bck.method==1
            [y2,f0s{ii},baseline(ii)]=calculatedDff0(Intensity(:,ii),1);
        elseif bck.method==2
            bck.fixPersti=[bck.grayNumber,qt(1)];
            [y2,f0s{ii},baseline(ii)]=calculatedDff0(Intensity(:,ii),bck.method,bck);
        elseif bck.method==3
            [y2,f0s{ii},baseline(ii)]=calculatedDff0(Intensity(:,ii),4);
        elseif bck.method==5
            bck.baseFixed=baseline_max(ii);
            [y2,f0s{ii},baseline(ii)]=calculatedDff0(Intensity(:,ii),bck.method,bck);
        end
        
        dff0(:,ii)=y2(:);
    end
    tt=toc;
    display(['geting intensities for individual ROIs takes ',num2str(tt,'%.1f'),'seconds'])
    file4=fileName(1).name;
    fileSave=['ROI_fj',strrep(fileName(1).name,'.tif','_Intensity.mat')];
    save(fullfile(filePath,fileSave),'Intensity', 'baseline', 'f0s', ...
        'dff0', 'stdImg','avg','bw','xy','p1','qt','order','file4','time');
    
    
else
    fileName2=dir(fullfile(filePath,['ROI_fj',fileName(1).name(id),'*.mat']));
    load(fullfile(filePath,fileName2(1).name));
end
