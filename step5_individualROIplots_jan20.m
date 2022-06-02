function resultPath=step5_individualROIplots(filePath,file,options) 
% plot results for individual ROIs
if nargin==0
%     filePath='E:\mydata_SC\2016_11_15\mouse354981\session01\stack';
    filePath='C:\Users\yajie.liang\Box\superior colliculus manuscript Liang\figures\Figure 1\fig1d_rawfiles\tuningmapandcurve';
%     filePath='Z:\Data\Yajie\20160928 342360 teststimu\s3_wz_ori\stack';
%     filePath='Z:\Data\Yajie\20160928 342360 teststimu\s5_rongw_ori_stripe\stack';
    file='ROI_fj*.mat';
    options.deleteN=[3,0];
    options.greenDot=1:options.deleteN(1);
    %% negative value issue
options.negativeMethod=0; % 1: set negative value as zero when calculating OSI,DSI,GOSI,GDSI; % 2: shift the curve; others: do nothing
    %% boundary
    
        options.LB=[1 -30 5 0.00001 5 -0.2]; %[amp1 theta sigma1 amp2 sigma2 dc];
        options.UB=[1.5 360 180 1.5 180 .5]; %[amp1 theta sigma1 amp2 sigma2 dc];
%         LB = [1 -30 5 0.00001 5 0]; % origonal value
%         UB = [1.5 360 180 1.5 180 .2];origonal value    

%         LB = [1 -30 5 0.00001 5 -0.2]; %Now allow the DC value to be negative for a better fitting. 
%         UB = [1.5 360 180 1.5 180 .5];
end
%%line432, fitting value is 70 or 15,20170510
%%line444, change fitting upper and lower bound,20170512
%%line1274, gDSI

options.file1=file;
options.trims=options.deleteN;

options.deleteFlag=0;%input.deleteROI is not working when input.deleteFlag is 0
options.deleteROI=[900:906];
  global minOptions
  minOptions=2;  %this is actually not used in this step
%         input.robustFit1.Flag=1;
%         input.robustFit1.dentrite=1;
output=main_ROI_process_version3_forGad2_protocol2_automated(filePath,options);
resultPath=output.resultPath{1};

% 
function input=main_ROI_process_version3_forGad2_protocol2_automated(filePath,input)

file1=input.file1;
fileName=dir(fullfile(filePath,file1));
%% processing stimulation; gad2Flag=0, turning curve; =1 normal; =2 no protocol;
p=length(fileName);
logFiles=cell(p,1);
id=1:11;
for ii=p:-1:1
    file1=fileName(ii).name;
    id2=id+regexp(file1,'\d\d\d\d\d\d\d\d_\d\d','once')-1;
    fileTmp=dir(fullfile(filePath,[file1(id2),'*.log']));
    p2=length(fileTmp);
    if p2>2
        fileName(ii)=[];
        logFiles(ii)=[];
    elseif p2==1
        logFiles{ii}=fullfile(filePath,fileTmp(1).name);
    end
end
p=length(logFiles);
% get log file info
f21_log=cell(p,1);
for ii=1:p
    if ~isempty(logFiles{ii})
        f21_log{ii} = f21log_load(logFiles{ii});
    end
    
end
[qt,gad2Flag]=obtainQuantity(f21_log);
p=length(fileName);
imgN=zeros(p,1);

for ii=1:p
         matObj=matfile(fullfile(filePath,fileName(ii).name));
        [m,n]=size(matObj,'Intensity');  
        imgN(ii)=m;
    if gad2Flag(ii)==1
    end
end
p=length(imgN);
repeatNO=imgN./qt(:,1)./qt(:,2);
repeatNO=floor(repeatNO);
[order,sequence]=getStimulusSequence(f21_log,repeatNO,gad2Flag);

%% deal with pupil

pupilPath=fullfile(filePath,'pupil');
if pupilPath==7
    folderFlag=1;
    
else
    [filePath0,name]=fileparts(filePath);
    pupilPath0=fullfile(filePath0,'pupil');
    if exist(pupilPath0)==7
        copyfile(pupilPath0,filePath)
        folderFlag=1;
    else
        folderFlag=0;
    
    end
    
end
pupilFlag=zeros(p,1);
eyeTrack=cell(p,1);
if folderFlag==1
    for ii=1:p
        file1=fileName(ii).name;
        id2=regexp(file1,'\d\d\d\d\d\d\d\d_\d\d');
        file2=file1(id2:id2+10);
%         fileName2=dir(fullfile(pupilPath,[file2,'*intensity_chosen.txt']));
        fileName2=dir(fullfile(pupilPath,[file2,'*area_chosen.txt']));
        if isempty(fileName2)
        else
            pupilFlag(ii)=1;
            if length(fileName2)>1
                fileName2=fileName2(1);

            end
            eyeTrack{ii}=load(fullfile(pupilPath,fileName2(1).name));            
        end

        
    end
end


p=length(fileName);
% for ii=p:-1:1
%     file1='comm_trash';
%     k1=strfind(fileName(ii).name,file1);
%     file1='dilate';
%     k2=strfind(fileName(ii).name,file1);
%     if ~isempty(k1) || ~isempty(k2)
%         fileName(ii)=[];    
%     end
% end
% 
% 
% % sectN=5:16;
% p=length(fileName);
fileAll=cell(p,1);
% timeStamp=regexp(file1,'\w\w\w\w-\w\w-\w\w_')
for ii=1:p
    ii
    trash=fileName(ii).name;
    sect1=regexp(trash,'20\w\w\w\w\w\w_\w');
    sectN=sect1+(0:11);
    sectString=fileName(ii).name(sectN);
%     trash=dir(fullfile(filePath,));
    fileAll{ii}=[sectString,'.mat'];
    
end

%% 
p=length(fileName);
ot_ds=cell(p,1);
result='ot_ds';
resultPath1=fullfile(filePath,result);
if exist(resultPath1)==7
    
else
    mkdir(resultPath1)
end

input.resultPath=cell(p,1);
for ii=1:p
    
    input.filePath=filePath;
    input.fileName=fileName(ii);
    if isempty(eyeTrack{ii})
        input.eyeTrack=eyeTrack{ii};
    else
        input.eyeTrack=eyeTrack{ii}(:,end);
    end
    
    input.pupilFlag=pupilFlag(ii);
%     input.eyeFlag=1;
    input.gad2Flag=gad2Flag(ii);
    input.order=order{ii};
    input.repeatNO=repeatNO(ii);
    input.test_sequence=sequence{ii};
    input.qt=[repeatNO(ii),qt(ii,2),qt(ii,1)];%22, 12,10
    
    file1=fileName(ii).name;
    id2=id+regexp(file1,'\d\d\d\d\d\d\d\d_\d\d','once')-1;
    fileTmp=dir(fullfile(filePath,[file1(id2),'*.mat']));

    
    if ~isempty(fileTmp)
            logFiles=fullfile(filePath,fileTmp(1).name);
        matObj=matfile(logFiles);
        if ~isempty(whos(matObj,'sequenceAngle'))
            orderTmp1=matObj.sequenceAngle;
            orderTmp2=orderTmp1(:);
            trial_angle=[matObj.trialNO,matObj.angleNO];
            repeatNO2=imgN(ii)./trial_angle(:,1)./trial_angle(:,2);
            repeatNO2=floor(repeatNO2);
            orderTmp3=ones(repeatNO2,1)*orderTmp2.';
            orderTmp4=orderTmp3(:);
            p2=length(orderTmp4);
              p1=p2;
                if p1<p2
                    orderTmp4=orderTmp4(1:p1);

                end 
                [Y,I]=sort(orderTmp4);
                order{ii}=[orderTmp4(:)*1,(1:length(I)).',Y(:)*1,I(:)]; 
                point_trial_angle=[repeatNO2,trial_angle];
                input.test_sequence=matObj.sequence+1;
                input.qt=point_trial_angle([1,3,2]);%22, 12,10
                input.order=order{ii};
                input.gad2Flag=0;
        end
        
    end   

    
    file=fileAll{ii};
    ROIfile=fileName(ii).name;
    [OT,DS,input,input.resultPath{ii}]=main_ROI_process_sub(file,filePath,ROIfile,input);
    ot_ds{ii}=[OT,DS];
%     file2=strrep(fileName(ii).name,'.tif','.mat');
%     save(fullfile(resultPath1,file2),'ot_ds')
    save(fullfile(resultPath1,[fileAll{ii},'_OTflag.txt']),'OT','-ascii');
    save(fullfile(resultPath1,[fileAll{ii},'_DSflag.txt']),'DS','-ascii');  
    if ii==p
        save(fullfile(resultPath1,'all.mat'),'ot_ds')
    end
%     =resultPath;
end




function [OT,DS,input,resultPath]=main_ROI_process_sub(file,filePath,ROIfile,input)

    input.ROIfileFlag=1;

input.ROIfile=ROIfile;

Datapath=filePath;
trims=input.trims;
projfns={file};

result='ROIs';
% result='trash';
resultPath1=fullfile(filePath,result);
if exist(resultPath1)==7
    
else
    mkdir(resultPath1)
end
result=[file(1:end-4)];
%     if input.ROIfileFlag==1
%         ROIfile=input.ROIfile;
%         result=[ROIfile(1:end-4)];
%     end
resultPath=fullfile(resultPath1,result);
if exist(resultPath)==7
    
else
    mkdir(resultPath)
end

for ifile=1:1:1   
    %Prepare experiment parameters ????
    orientation_pars=orientation_testpars(input);
    
    orientation_pars.Datapath=Datapath;
    orientation_pars.projfn=projfns{ifile};
    orientation_pars.ROIintensityfn=['ROI_' projfns{ifile}(1:end-4) '_Intensity.mat'];
    
    % 
    if input.ROIfileFlag==1
        orientation_pars.ROIintensityfn=ROIfile;
    end
    copyfile(fullfile(filePath,orientation_pars.ROIintensityfn),fullfile(resultPath,orientation_pars.ROIintensityfn))
    Datapath=resultPath;
    orientation_pars.Datapath=Datapath;
    savepath=Datapath;
    %
    orientation_pars.ROIsortdatafn=strrep(orientation_pars.ROIintensityfn,'_Intensity.mat','_Sortdata.mat');
    orientation_pars.ROItuningsfn=strrep(orientation_pars.ROIintensityfn,'_Intensity.mat','_Tunings_XXXX.mat'); %XXXX will be replaced by some data process parameters
    %%
    orientation_pars.testSeq=input.test_sequence;

    savefn=['0projInfo_' orientation_pars.projfn];
%     if ~exist(fullfile(savepath,savefn),'file')
        save(fullfile(savepath,savefn), 'orientation_pars', '-v7.3');
%     end
    
    %Load intensity data and sorting
    intensityData=load(fullfile(filePath,input.ROIfile));
    intensity=intensityData.Intensity;
    dff0=intensityData.dff0;
    f0s=intensityData.f0s;
    if isfield(input,'deleteFlag')
        if input.deleteFlag==1
            p=length(intensity);
            deleteROI=input.deleteROI;
            intensity(deleteROI,:)=nan;
            dff0(deleteROI,:)=nan;
            intensityData.Intensity(deleteROI,:)=nan;
            intensityData.dff0(deleteROI,:)=nan;
         
        end
    end    
    baseline=intensityData.baseline;
    
    [~,nROI]=size(intensityData.Intensity);
    qt=input.qt;
    sorted_intensity=reshape(intensityData.Intensity,qt(1),qt(2),qt(3),nROI);
    sorted_dff0=reshape(intensityData.dff0,qt(1),qt(2),qt(3),nROI);
    for irep=1:1:qt(3)
        ii=irep
        sorted_intensity(:,input.test_sequence(:,irep),irep,:)=sorted_intensity(:,:,irep,:);
        sorted_dff0(:,input.test_sequence(:,irep),irep,:)=sorted_dff0(:,:,irep,:);
    end
    input.intensity=intensity;
    input.dff0=dff0;
    input.f0s=f0s;
    input.sorted_intensity=sorted_intensity;
    input.sorted_dff0=sorted_dff0;
        save(fullfile(orientation_pars.Datapath,orientation_pars.ROIsortdatafn),...
            'intensity',...
            'dff0',...
            'f0s',...
            'baseline',...
            'sorted_intensity',...
            'sorted_dff0',...
            '-v7.3');
end

%%
%To intergrate frames for each test and generate normals
% Datapath='D:\LVI NTSR1 Orientation\';
Projinfo='0projInfo_*.mat';
Projinfofiles=dir(fullfile(Datapath, Projinfo));

%trims(1)+trims(2)+orientation_pars.pre_nframe_per_test<orientation_pars.nframe_per_test

for ifile=1:1:length(Projinfofiles)
    tuningfn=strrep(orientation_pars.ROItuningsfn,'Tunings_XXXX',sprintf('Tunings_TrimB%d_E%d',trims(1),trims(2)));
    if true(1,1)%~exist(fullfile(orientation_pars.Datapath,tuningfn),'file')

        [~,~,~,nROI]=size(sorted_dff0);
        
        sorted_data=sorted_dff0((orientation_pars.pre_nframe_per_test(1)+trims(1)+1):end-trims(2),:,:,:);
        data_mean_dir =squeeze(nanmean(sorted_data(:,:,:,:),1));
        DFFmean=squeeze(nanmean(data_mean_dir,2));
        DFFse=squeeze(nanstd(data_mean_dir,0,2)/sqrt(input.qt(3)));
        
        Pvalue_anova=nan(nROI,1);
        gOSI=NaN(nROI,1);
        gDSI=NaN(nROI,1);
        gDirection=NaN(nROI,1);
        
        for iROI=1:nROI
            [Pvalue_anova(iROI),~,~] = anova1((squeeze(nanmean(sorted_data(:,:,:,iROI),1)))',[], 'off');
            [gDSI(iROI),gDirection(iROI),~]=gDSI_cal_Kath(DFFmean(:,iROI), orientation_pars.testDatapoints(1),input);
            [gOSI(iROI),~]=gOSI_cal_kath(DFFmean(:,iROI),orientation_pars.testDatapoints(1),input);
        end
        
        save(fullfile(orientation_pars.Datapath,tuningfn),...
            'trims',...
            'DFFmean',...
            'DFFse',...
            'Pvalue_anova',...
            'gOSI',...
            'gDSI',...
            'gDirection',...
            '-v7.3');
        
        input.DFFmean=DFFmean;
        input.DFFse=DFFse;
        input.Pvalue_anova=Pvalue_anova;
        input.gOSI=gOSI;
        input.gDSI=gDSI;
        input.gDirection=gDirection;
        
    end
end
    
%%
% function ROI_Orientation_fitting(data_path)
% close all;
% clear all;
SSE_Threshold = 0.4;
RSquared_Threshold = 0.6;
data_path=resultPath;
% %layerVI CT NTSR1 neurons
% data_path = 'C:\data\Lu functional data\test2';
data_files = dir(fullfile(data_path, ['ROI_*_Tunings_TrimB',num2str(input.trims(1)),'_E',num2str(input.trims(2)),'.mat']));

nfiles = length(data_files);
% nCells = 0;

for ifile = 1:1:1
    
    fn=data_files(ifile).name;
    k1=strfind(fn,'ROI_');
    k2=strfind(fn, '_Tuning');
    projname_print=fn(k1+4:k2-1);
%     load(fullfile(data_path, data_files(ifile).name));
    [~,nROI] = size(DFFmean);
    
    OS_fitting_exitflag=-ones(nROI, 1);
    OS_fitting_sse=-ones(nROI, 1);
    OS_fitting_rsq=zeros(nROI, 1);
    
    OS_fitting_data_fitted_interp=[];
    OS_fitting_AMP1=NaN(nROI,1);
    OS_fitting_Theta=NaN(nROI,1);
    OS_fitting_Sigma=NaN(nROI,1);
    OS_fitting_AMP2=NaN(nROI,1);
    OS_fitting_DC=NaN(nROI,1);
    
    OSI = NaN(nROI,1);
    DSI=NaN(nROI,1);
    maxydata=NaN(nROI,1);
    
    %       for ii = 1:nROI
    %         nCells = nCells+1;
    global minOptions
    for ii = 1:nROI
        ydata=DFFmean(:,ii);
        SE=DFFse(:,ii);
        ydata = double(ydata);
        %%
        %deal with negative value: there is a replicate line at line 708
        if min(ydata)<0
            ydata=ydata-min(ydata);
            SE=SE-min(ydata);
            disp(['negative: ',num2str(ii)])
        end %remove negtive values
%         if min(ydata)<0
%             ydata(ydata<0)=0;
%             disp(['negative: ',num2str(ii)])
%         end %remove negtive values
        %%
        
        
        %normalize the data for fitting.
        maxydata(ii)=max(ydata);
        SE=SE/maxydata(ii);
        ydata=ydata/maxydata(ii);
        ydata=[ydata(end); ydata; ydata(1)];
        SE=[SE(end); SE; SE(1)];
        %to fit normal with a double guassian function
        %Initiate fitting parameters
        xdata=((0:numel(ydata)-1)-1)*(360/input.qt(2));
        xdata_interp=0:1:360;
        dc = min(ydata);
        [amp1 amp1_idx] = max(ydata);
        sigma1 = 70; %05/10/2017 Initial start value 15, then 70 for SC
        sigma2 = 70; %05/10/2017 Initial start value 15, then 70 for SC
        theta = xdata(amp1_idx);
        amp2_idx = amp1_idx-4;%amp2_idx = amp1_idx-6;
        if amp2_idx < 1
           amp2_idx = amp2_idx + 8;% amp2_idx = amp2_idx + 12;
        end
        amp2 = ydata(amp2_idx);
        
        start_point = [amp1 theta sigma1 amp2 sigma2 dc];
%         LB = [1 -30 5 0.00001 5 0];
%         UB = [1.5 360 180 1.5 180 .2];
%         LB = [1 -30 5 0.00001 5 -0.2]; %Now allow the DC value to be negative for a better fitting. 
%         UB = [1.5 360 180 1.5 180 .5];
        LB=input.LB;
        UB=input.UB;
        UB(6) = nanmean(ydata);
        
        %fitting
        [estimates, SSE, resnorm, residual, exitflag, lambda, jacobian] = fitcurve_OriGauss(xdata, ydata', start_point, LB, UB);
        
        %if first peak < second peak, exchange the theta's and redo
        %fitting ,which can make sue theta is always coressonding to
        %the biger peak at preferred orientation.
        if estimates(1)<estimates(4)
            tmpamp = estimates(1);
            estimates(1) = estimates(4);
            estimates(4) = tmpamp;
            estimates(2) = rem(estimates(2)+180, 360);
            [estimates, SSE, resnorm, residual, exitflag, lambda, jacobian] = fitcurve_OriGauss(xdata, ydata', estimates, LB, UB);
        end
        
        z = ori_gauss(estimates, xdata);
        z_interp = ori_gauss(estimates, xdata_interp);
        
        %Goodness of Fitting
        OS_fitting_exitflag(ii) = exitflag;
        %SSE
        OS_fitting_sse(ii) = SSE;
        %coefficient of determination, r_squred
        SSresid = sum(residual.^2);
        SStotal = (length(ydata)-1) * var(ydata);
        OS_fitting_rsq(ii) = 1 - SSresid/SStotal;
        
        %fittings
        data_fitted_interp(:,ii) = z_interp;
        if estimates(2)<0, estimates(2) = 360+estimates(2); end
        estimates(2) = rem(estimates(2),360);
        OS_fitting_AMP1(ii) = estimates(1);
        OS_fitting_Theta(ii) = estimates(2);
        OS_fitting_Sigma(ii) = estimates(3);
        OS_fitting_AMP2(ii) = estimates(4);
        OS_fitting_DC(ii) = estimates(6);
        
        z = ori_gauss(estimates, [estimates(2) rem(estimates(2)+180,360) rem(estimates(2)+90,360) rem(estimates(2)+270,360)]);
        OSI(ii) = (z(1)+z(2)-z(3)-z(4))/sum(z);
        DSI(ii)=(z(1)-z(2))/(z(1)+z(2));
        %             % % plot for visualization
        %             if SSE<=SSE_Threshold && OS_fitting_rsq(ii)>=RSquared_Threshold
        figure(2), clf
        ytheta=[0 OS_fitting_AMP1(ii)+OS_fitting_DC(ii)];
        xdirection=[gDirection(ii) gDirection(ii)]*180/pi;
        xtheta=[OS_fitting_Theta(ii) OS_fitting_Theta(ii)];
        hold on;
        errorbar(xdata(2:end-1), ydata(2:end-1), SE(2:end-1), 'ob');
        plot(xdata_interp, z_interp, '-c');
        plot(xtheta,ytheta,'r');
        plot(xdirection,ytheta,'g');
        hold off;
        xlim(gca, [xdata_interp(1) xdata_interp(end)]);
        set(gca, 'XTick', 0:(360/input.qt(2)):360);
        set(gca, 'TickDir', 'out');
        xlabel(gca, 'Degree(\circ)');
        ylabel(gca, 'Normalized \DeltaF/F');
        xlim(gca, [-5 365]);
        title(gca, ['Orentation normal, Pvalue=' num2str(Pvalue_anova(ii)) ', \theta=' num2str( estimates(2)) ', fiting SSE =' num2str(SSE) ', RSQ =' num2str(OS_fitting_rsq(ii)) ', OSI=' num2str(OSI(ii))]);
        yl=get(gca,'YLim');
        projname_print=strrep(projname_print,'_','-');
        text(5,yl(1)+0.05,[projname_print ':' num2str(ii)]);
        % % % % %                 print(gcf, fullfile(data_path,'OrientationFittingPlots.ps'), '-append');

    end
    
    OSfittingfn =strrep(data_files(ifile).name,'.mat','_OSfitting.mat');
    save(fullfile(data_path, OSfittingfn),...
        'maxydata',...
        'OS_fitting_AMP1',...
        'OS_fitting_Theta',...
        'OS_fitting_Sigma',...
        'OS_fitting_AMP2',...
        'OS_fitting_DC',...
        'OS_fitting_exitflag',...
        'OS_fitting_rsq',...
        'OS_fitting_sse',...
        'OSI',...
        'DSI',...
        '-v7.3');
end
% function [OT,DS]=ROI_Plot_Cells(data_path,input)

%combine all OSI for ROIs into single files
%%
% clear all;
% close all;

%Constants:
goodROI_threshold=5;
Pvalue_anova_threshold=0.05;
SSE_Threshold=0.4;
RSquared_Threshold=0.6;
FWHM_const=2.35482;

%L6 NTSR1
plot_celltype='L6 NTSR1';
% data_path = 'C:\data\Lu functional data\test2\';
% data_file='*Tunings_TrimB2_E0.mat';
data_file=['*Tunings_TrimB',num2str(input.trims(1)),'_E',num2str(input.trims(2)),'.mat'];
Tunings_files=dir(fullfile(data_path, data_file));

% intensity_files=dir([data_path strrep('*Tunings_TrimB2_E0.mat','Tunings_TrimB2_E0','_Intensity')]);
sortdata_files=dir(fullfile(data_path, strrep(data_file,data_file(2:end-4),'Sortdata')));
OSfitting_files=dir(fullfile(data_path, strrep(data_file,'.mat','_OSfitting.mat')));
% OSfitting_files=dir(fullfile((data_path, strrep('*Tunings_TrimB2_E0.mat','.mat','_OSfitting.mat'))));

nfiles = length(Tunings_files);

for ifile = 1:1:nfiles
    nROIs=0;
    matObj_sorted=matfile(fullfile(data_path, sortdata_files(ifile).name), 'Writable', true);
    matObj_tuning=matfile(fullfile(data_path, Tunings_files(ifile).name), 'Writable', true);
    matObj_fitting=matfile(fullfile(data_path, OSfitting_files(ifile).name), 'Writable', true);
    
    Intensity=matObj_sorted.intensity;    
    Baseline=matObj_sorted.baseline;
    DFF=matObj_sorted.dff0;
    F0s=matObj_sorted.f0s;
    Sorted_DFF=matObj_sorted.sorted_dff0;
       
    DFFmean =matObj_tuning.DFFmean;
    DFFse = matObj_tuning.DFFse;
    Pvalue_anova=matObj_tuning.Pvalue_anova;
    gOSI=matObj_tuning.gOSI;
    gDSI=matObj_tuning.gDSI;
    gDirection=matObj_tuning.gDirection;

    maxydata=matObj_fitting.maxydata;
    OS_fitting_AMP1=matObj_fitting.OS_fitting_AMP1;
    OS_fitting_AMP2=matObj_fitting.OS_fitting_AMP2;
    OS_fitting_DC=matObj_fitting.OS_fitting_DC;
    OS_fitting_Sigma=matObj_fitting.OS_fitting_Sigma;
    OS_fitting_Theta=matObj_fitting.OS_fitting_Theta;
    OS_fitting_exitflag=matObj_fitting.OS_fitting_exitflag;
    OS_fitting_rsq=matObj_fitting.OS_fitting_rsq;
    OS_fitting_sse=matObj_fitting.OS_fitting_sse;
    OSI=matObj_fitting.OSI;
    DSI=matObj_fitting.DSI;
    
    pdf_filename=strrep(OSfitting_files(ifile).name,'.mat','_rawplot.ps');
    if exist(fullfile(data_path,pdf_filename),'file')
        delete(fullfile(data_path,pdf_filename));
    end
    nROI=numel(OS_fitting_Theta);
    OT=zeros(nROI,1);
    DS=zeros(nROI,1);
    %plot each ROIs
    for iROI=1:1:nROI

        figure(1), clf;
        set(gcf,'PaperType','usletter','Unit','inches');
        set(gcf,'PaperPosition',[0.5 0.5 7.5 10]);
        ax(1)=subplot(20,4,1:12);
        plot_intensity(ax(1),Intensity(:,iROI),F0s{iROI});
        %% add static period
        qt=input.qt;
        qt23=qt(2)*qt(3);
        qt1=qt(1);
        yTrash=zeros(qt1,qt23);
        yTrash(input.greenDot,:)=1;
        yTrash(end-input.trims(1):end,:)=2;
         yTrash(input.trims(1)+1,:)=3;
        yTrash2=yTrash(:);
        xTrash=1:length(yTrash2);
        xTrash=xTrash(:);
        I_green=find(yTrash2==1);
        I_blue=find(yTrash2==2);
        I_red=find(yTrash2==3);
        set(ax(1),'NextPlot','add'); 
        plot(ax(1),xTrash(I_green),Intensity(I_green,iROI),'.g')
        plot(ax(1),xTrash(I_red),Intensity(I_red,iROI),'.r')
        xlabel(['green:',num2str(input.greenDot(1)),'-',num2str(input.greenDot(end)),';Red:',num2str(input.trims(1)+1)])
        
        %%
        if input.pupilFlag==1
           eyeTrack= input.eyeTrack;
           trash=Intensity(:,iROI);
           minI=min(trash);
           maxI=max(trash);
           eyeTrack2=(maxI-minI)/(1-0)*mat2gray(double(eyeTrack))+minI+(maxI-minI)*0.5;
           set(ax(1),'NextPlot','add');
            plot(ax(1),1:length(trash),eyeTrack2,'r')
            set(ax(1),'yLim',[minI,max(eyeTrack2(:))])
            

            
%                        eyeTrack= h.eyeTrack;
%            trash=Intensity(:,iROI);
%            minI=min(trash);
%            maxI=max(trash);
%            max_min=maxI-minI;
%            minI=minI-0.5*max_min;
%            maxI=maxI+0.5*max_min;
%            eyeTrack2=(maxI-minI)/(1-0)*mat2gray(double(eyeTrack))+minI+(maxI-minI)*0;
%            set(ax(1),'NextPlot','add');
%             h_trash=plot(ax(1),1:length(trash),eyeTrack2,'r')
%             set(ax(1),'yLim',[minI,maxI])
%             uistack(h_trash,'bottom') 
            
%             legend()
        end
        set(ax(1),'NextPlot','replace'); 
        if iROI==6
            trash=0;
        end
        ax(2)=subplot(20,4,17:28);
        plot_dff(ax(2),DFF(:,iROI),F0s{iROI});
        txt_title=sortdata_files(ifile).name(5:end);
        txt_title=strrep(txt_title,'_Sortdata.mat',':');
        txt_title=strrep(txt_title,'_','-');
        title(ax(1),[txt_title 'ROI#' sprintf('%03d', iROI)]);
        %%
         set(ax(2),'NextPlot','add'); 
        plot(ax(2),xTrash(I_red),DFF(I_red,iROI),'.r')
        xlabel(['Red:',num2str(input.trims(1)+1)])
         set(ax(1),'NextPlot','replace'); 
%         plot(ax(2),xTrash(I_blue),Intensity(I_blue,iROI),'.b')
        
        ax(3)=subplot(20,4,[37 38 41 42 45 46 49 50 53 54]);
        plot_tuningcurve(ax(3),Sorted_DFF(:,:,:,iROI),input.qt,input.trims,'normal',10);%%trims was [5 0]
        %%
%         plot_tuningcurve(ax(3),Sorted_DFF(:,:,:,iROI),[29 8 8],[5 0],'normal',10);
%         plot_tuningcurve(ax(3),Sorted_DFF(:,:,:,iROI),[22 12 10],[5 0],'normal',10);
        
        ax(4)=subplot(20,4,[39 40 43 44 47 48 51 52 55 56]);
        estimates=[OS_fitting_AMP1(iROI) OS_fitting_Theta(iROI) OS_fitting_Sigma(iROI) OS_fitting_AMP2(iROI) OS_fitting_Sigma(iROI) OS_fitting_DC(iROI)];
        xdata_interp=0:1:360;
        ydatafit=ori_gauss(estimates, xdata_interp);
        ytheta=[0 OS_fitting_AMP1(iROI)+OS_fitting_DC(iROI)];
        xdirection=[gDirection(iROI) gDirection(iROI)]*180/pi;
        xtheta=[OS_fitting_Theta(iROI) OS_fitting_Theta(iROI)];
        xdata=(0:(input.qt(2)-1))*360/(input.qt(2));
        ydata=DFFmean(:,iROI);
        ydatase=DFFse(:,iROI);
        
        if min(ydata)<0 %this is a replicate of line 430
            ydata=ydata-min(ydata);
            ydatase=ydatase-min(ydata);
        end
        
%         if min(ydata)<0
%             ydata(ydata<0)=0;
%             disp(['negative: ',num2str(ii)])
%         end %remove negtive values
        
        hold on;
        errorbar(xdata, ydata/maxydata(iROI), ydatase/maxydata(iROI), 'ob','Parent',ax(4));
        plot(ax(4),xdata_interp, ydatafit, '-c');
        plot(ax(4),xtheta,ytheta,'r');
        plot(ax(4),xdirection,ytheta,'--g');
        hold off;
        xlim(ax(4), [xdata_interp(1) xdata_interp(end)]);
        set(ax(4), 'XTick', 0:360/(input.qt(2)):360);
        set(ax(4), 'TickDir', 'out');
        xlabel(ax(4), 'Degree(\circ)');
        ylabel(ax(4), 'Normalized \DeltaF/F');
        xlim(ax(4), [-5 365]);
        title(ax(4), ['Pvalue=' num2str(Pvalue_anova(iROI)) ', \theta=' num2str( estimates(2)) ', SSE =' num2str(OS_fitting_sse(iROI)) ', RSQ =' num2str(OS_fitting_rsq(iROI))]);
        yl=get(ax(4),'YLim');
        
         ax(5)=subplot(20,4,[63 64 67 68 71 72 75 76 79 80]);
         dori = diff(xdata(1:2));
         xdata = rem(xdata+90,360)/180*pi;
         mw = max(max(ydata/maxydata(iROI)), max(ydatafit));
         r = circ_r(xdata,(ydata/maxydata(iROI))',dori) * mw;
         phi = circ_mean(xdata,(ydata/maxydata(iROI))');
         %     xdata_interp = xdata_interp/180*pi;
         hold on
         zm = r*exp(i*phi');
         
         % draw a unit circle
         zz = exp(i*linspace(0, 2*pi, 101)) * mw;
         plot(real(zz),imag(zz),'-', 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5])
         plot(real(zz/2),imag(zz/2),'-','LineWidth', 0.5, 'Color', [0.5 0.5 0.5])
         plot(real(0/2),imag(0/2),'-','LineWidth', 0.5, 'Color', [0.5 0.5 0.5])
         plot([-mw mw], [0 0], '-', [0 0], [-mw mw], '-', 'Color', [0.5 0.5 0.5])
         %
         polar([xdata xdata(1)], [(ydata/maxydata(iROI))' ydata(1)/maxydata(iROI)], 'ok')
         polar([xdata_interp+90 xdata_interp(1)+90]/180*pi, [ydatafit ydatafit(1)], 'b');
         hold off
         
         formatSubplot(gca,'ax','square','box','off','lim',[-mw-0.1 mw+0.1 -mw-0.1 mw+0.1]);
         set(gca,'xtick',[]);
         set(gca,'ytick',[]);
         text(-max(real(zz))/20,max(real(zz))*1.1,'Up');
         text(-max(real(zz))/9,-max(real(zz)*1.1),'Down');
         text(max(real(zz))*1.1,max(real(zz))/5,'Posterior','Rotation',270);
         text(-max(real(zz))*1.1,-max(real(zz))/6,'Anterior','Rotation',90);
         axis off;
         
         ax(6)=subplot(20,4,[61 62 65 66 69 70 73 74 77 78]);
         set(ax(6),'XLim',[0 1],'Ylim',[0 1]);
         axis off;
         
         is_goodROI='No';
         if (maxydata(iROI)>=goodROI_threshold), is_goodROI='Yes'; end
         text(0,0.95,sprintf('This ROI is a good ROI(Max respons>=%d%%): %s\n max response=%.2f',goodROI_threshold,is_goodROI,maxydata(iROI)),'Parent',ax(6));
         isgoodfitting='No';
         if (OS_fitting_sse(iROI)<=SSE_Threshold)&&...
            (OS_fitting_rsq(iROI)>=RSquared_Threshold)&&...
            (maxydata(iROI)>=goodROI_threshold)
            isgoodfitting='Yes';
         end
         text(0,0.8,sprintf('A good fitting(fitting SSE<%.2f AND RSQ>%.2f): %s',SSE_Threshold,RSquared_Threshold,isgoodfitting),'Parent',ax(6));
         
         isOrientationTuned='No';
         if (Pvalue_anova(iROI)<Pvalue_anova_threshold)&&strcmp(isgoodfitting,'Yes'),isOrientationTuned='Yes'; end
%          if (Pvalue_anova(iROI)<Pvalue_anova_threshold); end %20210113
         text(0,0.6,sprintf('This ROI is OS(anova test pValue<%.3f): %s',Pvalue_anova_threshold,isOrientationTuned),'Parent',ax(6));
         osstr=sprintf('OS: Theta=%.2f, FWHM=%.2f,\n gOSI=%.2f, OSI=%.2f, DSI=%.2f',estimates(2),estimates(3)*FWHM_const,gOSI(iROI),OSI(iROI),DSI(iROI));
         if strcmp(isOrientationTuned,'No')
%            osstr=sprintf('OS: Theta=%.2f, FWHM=%.2f,\n gOSI=%.2f, OSI=%.2f, DSI=%.2f',NaN,NaN,NaN,NaN,NaN);
           osstr=sprintf('OS: Theta=%.2f, FWHM=%.2f,\n gOSI=%.2f, OSI=%.2f, DSI=%.2f',estimates(2),estimates(3)*FWHM_const,gOSI(iROI),OSI(iROI),DSI(iROI));
         else
             OT(iROI)=1;
         end
         text(0,0.45,osstr,'Parent',ax(6));
         
         isDS='No';
         if (DSI(iROI)>=0.5)&&strcmp(isOrientationTuned,'Yes'), isDS='Yes'; end
         text(0,0.3,sprintf('This ROI is DS(DSI>=0.5): %s',isDS),'Parent',ax(6));
         DSstr=sprintf('DS: Theta=%.2f, gDirection=%.2f, FWHM=%.2f,\n gDSI=%.2f, DSI=%.2f',estimates(2),gDirection(iROI)*180/pi,estimates(3)*FWHM_const,gDSI(iROI),DSI(iROI));
         if strcmp(isDS,'No')
%              DSstr=sprintf('DS: Theta=%.2f, gDirection=%.2f, FWHM=%.2f,\n gDSI=%.2f, DSI=%.2f',NaN,NaN,NaN,NaN,NaN);
             DSstr=sprintf('DS: Theta=%.2f, gDirection=%.2f, FWHM=%.2f,\n gDSI=%.2f, DSI=%.2f',estimates(2),gDirection(iROI)*180/pi,estimates(3)*FWHM_const,gDSI(iROI),DSI(iROI));
         else
            DS(iROI)=1;
         end
         text(0,0.15,DSstr,'Parent',ax(6));
         
        
        %print this figure into a PDF file
        %print this figure into a PDF file
        if iROI==1
            if exist(fullfile(data_path,pdf_filename))==2
                delete(fullfile(data_path,pdf_filename))
            end
         
        end
        print(gcf,'-dpsc2',fullfile(data_path,pdf_filename),'-append');
%         savefig(strrep(gcf,'.eps','.fig'));
        saveas(gcf,fullfile(data_path,[pdf_filename(4+(1:11)),'_',num2str(iROI,'%03d'),'.eps'],''))
        saveas(gcf,fullfile(data_path,[pdf_filename(4+(1:11)),'_',num2str(iROI,'%03d'),'.eps','.fig'],''))
%         
        filePath='C:\Users\lur\Dropbox (HHMI)\psfiles';
        if exist(filePath)==7

        if iROI==1
            if exist(fullfile(filePath,pdf_filename))==2
                delete(fullfile(filePath,pdf_filename))
            end
         
        end 
            print(gcf,'-dpsc2',fullfile(filePath,pdf_filename),'-append');        
        end
    end
    file1=strrep(matObj_sorted.Properties.Source,'.mat','_OT.txt');
    file2=strrep(matObj_sorted.Properties.Source,'.mat','_DS.txt');
    save(file1,'OT','-ascii');
    save(file2,'DS','-ascii');
    
    
    
end
% step6([data_path,'\'])
function output = plot_tuningcurve(ahandle,sorted_data,ndimens,trims,plot_title,nblanks)

if nargin<6, nblanks=10;end
if nargin<5, plot_title='normal'; end
if nargin<4, trims=[0 0]; end

output = 0;

data_plot = nan(ndimens(1)+nblanks, ndimens(2), ndimens(3));
data_plot(1:ndimens(1),:,:) = sorted_data;
data_mean_rep = reshape(data_plot, ndimens(2)*(ndimens(1)+nblanks), ndimens(3));
data_se_rep = nanstd(data_mean_rep,0,2)/sqrt(ndimens(3));
data_mean_rep = nanmean(data_mean_rep,2);

data_mean_dir =squeeze(nanmean(sorted_data(trims(1)+1:end-trims(2),:,:),1));
ydata_mean=nanmean(data_mean_dir,2);
yerror_mean=nanstd(data_mean_dir,0,2)/sqrt(ndimens(3));
if min(ydata_mean)<0
    data_mean_dir=data_mean_dir-min(yerror_mean);
    yerror_mean=yerror_mean-min(yerror_mean);
    yerror_mean=yerror_mean-min(yerror_mean);
end

hold on
xdata = 1:1:ndimens(2)*(ndimens(1)+nblanks);
errorbar(xdata(1:end-nblanks+1), data_mean_rep(1:end-nblanks+1), data_se_rep(1:end-nblanks+1),'Color',[0.7 0.7 0.7],'LineWidth',0.5, 'Parent', ahandle);

xdata_mean=(0:ndimens(2)-1)*(ndimens(1)+nblanks)+round(ndimens(1)/2);
errorbar(ahandle,xdata_mean,ydata_mean,yerror_mean,'LineStyle','none','Marker', 'o','MarkerSize',10,'LineWidth',2, 'Parent', ahandle);
hold off;

box('off');
set(ahandle, 'TickDir', 'out');
set(ahandle, 'XLim', [xdata(1)-25 xdata(end)+25],'XTick', xdata_mean, 'XTickLabel',(1:ndimens(2))-1);
set(ahandle, 'FontName', 'arial', 'FontSize', 10);

ylabel(ahandle,'\DeltaF/F0','FontName','arial','FontAngle','italic','FontSize',12);
xlabel(ahandle,'Test #','FontName','arial','FontSize',12);
function z = ori_gauss(params, xdata)
%use this double gaussian function to fit orientation tuning curve
%     amp1= params(1);
%     theta = params(2);
%     sigma1 = params(3);
%     amp2 = params(4);
%     sigma2 = params(5);
%     dc = params(6);
% created by Wenzhi, 04/10/2013

%different sigma
% z =params(1)*exp(-((-360*(sign(xdata-params(2)-180)>0)+360*(sign(xdata-params(2)+180)<0)+xdata-params(2)).^2)/params(3)/params(3)/2)...
%     +params(4)*exp(-((360*(sign(xdata-params(2))<0)+xdata-params(2)-180).^2)/params(5)/params(5)/2)...
%     +params(6);

%same sigma
z =params(1)*exp(-((-360*(sign(xdata-params(2)-180)>0)+360*(sign(xdata-params(2)+180)<0)+xdata-params(2)).^2)/params(3)/params(3)/2)...
    +params(4)*exp(-((360*(sign(xdata-params(2))<0)+xdata-params(2)-180).^2)/params(3)/params(3)/2)...
    +params(6);
function [estimates, SSE, resnorm, residual,exitflag, lambda, jacobian] = fitcurve_OriGauss(xdata, ydata, start_point, LB, UB)
opt = optimset('Display','off');

% opt=optimset('MaxFunEvals', 90000000, 'MaxIter', 900000000, 'TolFun',  1.000000e-08, 'Display','off');
[estimates,resnorm, residual,exitflag, output, lambda, jacobian] = lsqcurvefit(@ori_gauss, start_point, xdata, ydata, LB, UB, opt);
z = ori_gauss(estimates, xdata);
SSE = sum((z-ydata).*(z-ydata));
function [output_args]=plot_intensity(ahandle,intensity,f0s)
%Inputs:    ahandle
%           intensity and indics of data points selected as F0 as baseline
%           of dff calculation.
%created by Wenzhi Sun, 08.17.2015
if nargin<3, f0s=[]; end
if ishandle(ahandle)
    cla(ahandle);
    Imax = max(intensity);
    Imin = min(intensity);
    xdata = 1:length(intensity);
    hold on
    plot(ahandle,xdata, intensity, 'k', 'LineWidth', .5);
%     plot(ahandle,f0s, intensity(f0s), 'm.', 'MarkerSize', 6, 'MarkerFaceColor', 'none');
    hold off;
    box('off');
    set(ahandle, 'TickDir', 'out');
    set(ahandle, 'XLim', [xdata(1)-25 xdata(end)+25]);
    set(ahandle, 'YLim', [Imin*0.95 Imax+Imin/20]);
    set(ahandle, 'FontName', 'arial', 'FontSize', 10);
    
    ylabel(ahandle,'Fluo. Intensity (a.u.)','FontName','arial','FontAngle','normal','FontSize',12);
    xlabel(ahandle,'Frame #','FontName','arial','FontSize',12);
end
function [output_args]=plot_dff(ahandle,dff,f0s)
%Inputs:    ahandle
%           dff and indics of f0s 
%created by Wenzhi Sun, 08.17.2015
if nargin<3, f0s=[]; end
if ishandle(ahandle)
    cla(ahandle);
    dmax = max(dff);
    dmin = min(dff);
    xdata = 1:length(dff);
    hold on
    plot(ahandle,xdata, dff, 'k', 'LineWidth', .5);
    plot(ahandle,f0s, dff(f0s), 'bx', 'MarkerSize', 0.5, 'MarkerFaceColor', 'none');
    hold off;
    box('off');
    set(ahandle, 'TickDir', 'out');
    set(ahandle, 'XLim', [xdata(1)-25 xdata(end)+25]);
    set(ahandle, 'YLim', [dmin*0.95 dmax+dmin/20]);
    set(ahandle, 'FontName', 'arial', 'FontSize', 10);
    
    ylabel(ahandle,'\DeltaF/F0','FontName','arial','FontAngle','italic','FontSize',12);
    xlabel(ahandle,'Frame #','FontName','arial','FontSize',12);  
end

function r = circ_r(alpha, w, d, dim)
% r = circ_r(alpha, w, d)
%   Computes mean resultant vector length for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [w		number of incidences in case of binned angle data]
%     [d    spacing of bin centers for binned data, if supplied 
%           correction factor is used to correct for bias in 
%           estimation of r, in radians (!)]
%     [dim  compute along this dimension, default is 1]
%
%     If dim argument is specified, all other optional arguments can be
%     left empty: circ_r(alpha, [], [], dim)
%
%   Output:
%     r		mean resultant length
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N.I. Fisher
%   Topics in circular statistics, S.R. Jammalamadaka et al. 
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if nargin < 4
  dim = 1;
end

if nargin < 2 || isempty(w) 
  % if no specific weighting has been specified
  % assume no binning has taken place
	w = ones(size(alpha));
else
  if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1) 
    error('Input dimensions do not match');
  end 
end

if nargin < 3 || isempty(d)
  % per default do not apply correct for binned data
  d = 0;
end

% compute weighted sum of cos and sin of angles
r = sum(w.*exp(1i*alpha),dim);

% obtain length 
r = abs(r)./sum(w,dim);

% for data with known spacing, apply correction factor to correct for bias
% in the estimation of r (see Zar, p. 601, equ. 26.16)
if d ~= 0
  c = d/2/sin(d/2);
  r = c*r;
end

function [mu ul ll] = circ_mean(alpha, w, dim)
%
% mu = circ_mean(alpha, w)
%   Computes the mean direction for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [w		weightings in case of binned angle data]
%     [dim  compute along this dimension, default is 1]
%
%     If dim argument is specified, all other optional arguments can be
%     left empty: circ_mean(alpha, [], dim)
%
%   Output:
%     mu		mean direction
%     ul    upper 95% confidence limit
%     ll    lower 95% confidence limit 
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N. I. Fisher
%   Topics in circular statistics, S. R. Jammalamadaka et al. 
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if nargin < 3
  dim = 1;
end

if nargin < 2 || isempty(w)
  % if no specific weighting has been specified
  % assume no binning has taken place
	w = ones(size(alpha));
else
  if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1) 
    error('Input dimensions do not match');
  end 
end

% compute weighted sum of cos and sin of angles
r = sum(w.*exp(1i*alpha),dim);

% obtain mean by
mu = angle(r);

% confidence limits if desired
if nargout > 1
  t = circ_confmean(alpha,0.05,w,[],dim);
  ul = mu + t;
  ll = mu - t;
end

function t = circ_confmean(alpha, xi, w, d, dim)
%
% t = circ_mean(alpha, xi, w, d, dim)
%   Computes the confidence limits on the mean for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [xi   (1-xi)-confidence limits are computed, default 0.05]
%     [w		number of incidences in case of binned angle data]
%     [d    spacing of bin centers for binned data, if supplied 
%           correction factor is used to correct for bias in 
%           estimation of r, in radians (!)]
%     [dim  compute along this dimension, default is 1]
%
%   Output:
%     t     mean +- d yields upper/lower (1-xi)% confidence limit
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N. I. Fisher
%   Topics in circular statistics, S. R. Jammalamadaka et al. 
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if nargin < 5
  dim = 1;
end

if nargin < 4 || isempty(d)
  % per default do not apply correct for binned data
  d = 0;
end

if nargin < 3 || isempty(w)
  % if no specific weighting has been specified
  % assume no binning has taken place
	w = ones(size(alpha));
else
  if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1) 
    error('Input dimensions do not match');
  end 
end

% set confidence limit size to default
if nargin < 2 || isempty(xi)
  xi = 0.05;
end

% compute ingredients for conf. lim.
r = circ_r(alpha,w,d,dim);
n = sum(w,dim);
R = n.*r;
c2 = chi2inv((1-xi),1);

% check for resultant vector length and select appropriate formula
t = zeros(size(r));

for i = 1:numel(r)
  if r(i) < .9 && r(i) > sqrt(c2/2/n(i))
    t(i) = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));  % equ. 26.24
  elseif r(i) >= .9
    t(i) = sqrt(n(i)^2-(n(i)^2-R(i)^2)*exp(c2/n(i)));      % equ. 26.25
  else 
    t(i) = NaN;
    warning('Requirements for confidence levels not met.');
  end
end

% apply final transform
t = acos(t./R);
function formatSubplot(handle,varargin)

args.fs = 7;
args.xl = [];
args.yl = [];
args.box = [];
args.ax = [];
args.lim = [];
args.tt = [];
args.xt = [];
args.yt =[];
args = parseVarArgs(args,varargin{:});

set(handle,'fontsize',args.fs)
if ~isempty(args.xl)
  xlabel(handle,args.xl)
end
if ~isempty(args.yl)
  ylabel(handle,args.yl)
end
if ~isempty(args.tt)
  title(handle,args.tt)
end
if ~isempty(args.box)
  set(handle,'box',args.box)
end
if ~isempty(args.ax)
  axis(handle,args.ax)
end
if ~isempty(args.lim)
  axis(handle,args.lim)
end
if ~isempty(args.yt)
  set(handle,'ytick',args.yt)
end
if ~isempty(args.xt)
  set(handle,'xtick',args.xt)
end
function params = parseVarArgs(params,varargin)
% Parse variable input arguments supplied in name/value format.
%
%    params = parseVarArgs(params,'property1',value1,'property2',value2) sets
%    the fields propertyX in p to valueX.
%
%    params = parseVarArgs(params,varargin{:},'strict') only sets the field
%    names already present in params. All others are ignored.
%
% AE 2007-06-01

if isempty(varargin)
    return
end

% check if correct number of inputs
if mod(length(varargin),2)
    if ~strcmp(varargin{end},'strict')
        err.message = 'Name and value input arguments must come in pairs.';
        err.identifier = 'parseVarArgs:wrongInputFormat';
        error(err)
    else
        % in 'strict' case, remove all fields that are not already in params
        fields = fieldnames(params);
        ndx = find(~ismember(varargin(1:2:end-1),fields));
        varargin([2*ndx-1 2*ndx end]) = [];
    end
end

% parse arguments
for i = 1:2:length(varargin)
    if ischar(varargin{i})
        params.(varargin{i}) = varargin{i+1};
    else
        err.message = 'Name and value input arguments must come in pairs.';
        err.identifier = 'parseVarArgs:wrongInputFormat';
        error(err)
    end
end

function [order,sequence]=getStimulusSequence(f21_log,repeatNO,gad2Flag)

p=length(f21_log);

order=cell(p,1);
sequence=cell(p,1);
for ii=1:p
    ii
%     ii
    if gad2Flag(ii)==1 || gad2Flag(ii)==2
        sequence{ii}=[0,3,4,6,6,0,4,6;5,4,2,1,7,2,6,4;3,7,7,7,1,4,1,1;4,5,1,5,3,6,7,2;7,1,3,2,5,1,5,7;1,6,5,3,4,3,3,0;2,2,6,0,2,7,2,5;6,0,0,4,0,5,0,3];
        orderTmp1=45*sequence{ii};
        orderTmp2=orderTmp1(:);  
        sequence{ii}=sequence{ii}+1;
    elseif gad2Flag(ii)==0
        orderTmp2=f21_log{ii}.test_value_sequence(:);
        sequence{ii}=f21_log{ii}.test_sequence;
        sequence{ii}=sequence{ii}+1;

    end
    
    orderTmp3=ones(repeatNO(ii),1)*orderTmp2.';
    orderTmp4=orderTmp3(:);
    % for future use
          p2=length(orderTmp4);
          p1=p2;
            if p1<p2
                orderTmp4=orderTmp4(1:p1);
            end 
            [Y,I]=sort(orderTmp4);
            order{ii}=[orderTmp4(:)*1,(1:length(I)).',Y(:)*1,I(:)];        
end
function [qt,gad2Flag]=obtainQuantity(f21_log)

p2=length(f21_log);
qt=zeros(p2,2);
gad2Flag=zeros(p2,1);
for ii=1:p2
    if ~isempty(f21_log{ii})
        if f21_log{ii}.isvalid==1
            qt(ii,:)=[f21_log{ii}.repeats,f21_log{ii}.data_points];
        else
            qt(ii,:)=[8,8];
            gad2Flag(ii)=1;

        end
    else
        qt(ii,:)=[8,8];gad2Flag(ii)=1;
    end    

end


function orientation_pars_template=orientation_testpars(input)
    
    orientation_pars_template.testType='normal';
    orientation_pars_template.tc_mode='Ori_CRF_Standard';
    if input.gad2Flag==1
        orientation_pars_template.nframe_per_test=input.qt(1);
        orientation_pars_template.pre_nframe_per_test=[10];
        orientation_pars_template.post_nframe_per_test=[0];
        orientation_pars_template.testSeq=[];
        orientation_pars_template.testDatapoints=input.qt(2);
        orientation_pars_template.testReps=input.qt(3);         
    else
    orientation_pars_template.nframe_per_test=input.qt(1);
    orientation_pars_template.pre_nframe_per_test=[0];
    orientation_pars_template.post_nframe_per_test=[0];
    orientation_pars_template.testSeq=[];
    orientation_pars_template.testDatapoints=input.qt(2);
    orientation_pars_template.testReps=input.qt(3);     
    end

    
    orientation_pars_template.logfn='';
    orientation_pars_template.projfn='';
    orientation_pars_template.ROIintensityfn='';
    orientation_pars_template.ROIsortdatafn='';
    orientation_pars_template.ROItuningsfn='';
    orientation_pars_template.Datapath='';
    orientation_pars_template.iferror=0;
    orientation_pars_template.errorType='';
function [gDSI,pref_direction,exception]=gDSI_cal(ydata,datapoints)
% A direction-selectivity index (DSI) as the vector sum of responses
% normalized by the scalar sum of responses (such that the index varies
% between 0 and 1).
% The angle of this vector sum defined the preferred direction of each
% cell.
%created by Wenzhi Sun, July 27 2015
if isvector(ydata)
    if eq(numel(ydata),datapoints)
        if iscolumn(ydata), ydata=ydata'; end
        theta=0:(2*pi/datapoints):(2*pi-0.1);
        %for 8 datapoints-measurement theta=[0 pi/4 pi/2 pi3/4 pi pi5/4 pi3/2 and pi7/4]
%         if min(ydata)<0, ydata=ydata-min(ydata); end %remove the negtive value

        gDSI=abs(sum(ydata.*exp(1i.*theta))/sum(ydata));
        pref_direction=angle(sum(ydata.*exp(1i.*theta)));
        %P = angle(Z) returns the phase angles, in radians, for each element of complex array Z. The angles lie between ±pi.
        %to wrap the phase angle to [0-2pi)
        if pref_direction<0, pref_direction=pref_direction+pi*2; end
        exception=[];
    else
        gDSI=nan;
        pref_direction=nan;
        exception=MException(['ws_IMAGEBOX::' 'THE DATA POINTS OF INPUT DATA DOSN''T MATCH THE INPUT DATAPOINTS!'],'Iputs Not Match');
    end
else
    gDSI=nan;
    pref_direction=nan;
    exception=MException(['ws_IMAGEBOX::' 'THE FIRST INPUT TO FUNCTION gDSI_cal MUST BE A VECTOR!'],'Wrong input type');   
end    
% return;
function [gOSI,exception]=gOSI_cal(ydata,datapoints)
%Inputs: ydata a vector of responseses and how many datapoints in this
%measurement. numel(ydata)==datapoints
%created by Wenzhi Sun, July 27 2015
    if isvector(ydata)
        if eq(numel(ydata),datapoints)
            if iscolumn(ydata), ydata=ydata'; end
            theta=0:(2*pi/datapoints):(2*pi-0.1);
%             if min(ydata)<0, ydata=ydata-min(ydata); end %remove the negtive value
            
            gOSI=abs(sum(ydata.*exp(1i*2*theta))/sum(ydata));
            exception=[];
        else
            gOSI=nan;
            exception=MException(['ws_IMAGEBOX::' 'THE DATA POINTS OF INPUT DATA DOSN''T MATCH THE INPUT DATAPOINTS!'],'Iputs Not Match');
        end
    else
        gOSI=nan;
        exception=MException(['ws_IMAGEBOX::' 'THE FIRST INPUT TO FUNCTION gOSI_cal MUST BE A VECTOR!'],'Wrong input type');
    end
% end




function f21_log = f21log_load(logfile)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Load and prorcess the log file
if ~isequal(exist(logfile), 2)
    disp('Can not find the log file, the program will stop!!');
    return;
end

lines = f21log_get(logfile);
if ~iscellstr(lines)
   % errordlg('Can''t load log file data correctly.','File Load Error','replace')
   f21_log.isvalid = 0;
   return;
end

% search log file data
file_type = f21log_search(lines,'FileInfo','FileType','JQK');
if ~strcmp(file_type,'f21lv test log file')
%    % errordlg('This is not a valid f21 log file.','File Type Error','replace')
   f21_log.isvalid = 0;
   return;
end

f21_log.file_type = file_type;

file_version = f21log_search(lines,'FileInfo','FileVersion','JQK');
if ~strcmp(file_version,'2.0')
%    % errordlg('Only version 2.0 log file is supported.','File Version Error','replace')
   f21_log.isvalid = 0;
   return;
else
    f21_log.isvalid = 1;
end
f21_log.file_version = file_version;

file_time = f21log_search(lines,'FileInfo','FileTime','NULL');
f21_log.stop_time = file_time;

test_type = f21log_search(lines,'TestInfo','TestType','JQK');
if ~strcmp(test_type,'tuning curve')
   % % errordlg('This is not a tuning curve test log file.','Test Type Error','replace')
   f21_log.isvalid = 0;
   return;
end
f21_log.test_type = test_type;

test_name = f21log_search(lines,'TestInfo','Testname','JQK');
f21_log.test_name = test_name;

str_temp = f21log_search(lines,'TestInfo','RemoteRefreshRate','1');
remote_refresh_rate = sscanf(str_temp,'%f');
f21_log.remote_refresh_rate = remote_refresh_rate;

str_temp = f21log_search(lines,'TestInfo','RefreshRate','1');
refresh_rate = sscanf(str_temp,'%f');
f21_log.refresh_rate = refresh_rate;

% use 0.1 ms as an unit comparing to 10khz acq rate;
UNIT_PER_SEC = 1000*10;
unit_per_frame = UNIT_PER_SEC/refresh_rate;
f21_log.unit_per_frame = unit_per_frame;

tc_mode = f21log_search(lines,'TunningCurve','TCMode','JQK');
if strcmp(tc_mode,'JQK')
   % errordlg('This is not a new-style tuning curve test log file.','TC Type Error','replace')
   f21_log.isvalid = 0;
   return
end
f21_log.tc_mode = tc_mode;

gratings_type = f21log_search(lines,'TunningCurve','CRFShowType','JQK');
f21_log.gratings_type = gratings_type;

str_temp = f21log_search(lines,'TunningCurve','OriX','9999');
ori_x = sscanf(str_temp,'%f');
f21_log.ori_x = ori_x;

str_temp = f21log_search(lines,'TunningCurve','OriY','9999');
ori_y = sscanf(str_temp,'%f');
f21_log.ori_y = ori_y;

str_temp = f21log_search(lines,'TunningCurve','GratingsW','9999');
gratings_w = sscanf(str_temp,'%f');
f21_log.gratings_w = gratings_w;

str_temp = f21log_search(lines,'TunningCurve','GratingsH','9999');
gratings_h = sscanf(str_temp,'%f');
f21_log.gratings_h = gratings_h;

str_temp = f21log_search(lines,'TunningCurve','GratingsSpatialFrequency','9999');
gratings_spatial_frequency = sscanf(str_temp,'%f');
f21_log.gratings_spatial_frequency = gratings_spatial_frequency;

str_temp = f21log_search(lines,'TunningCurve','GratingsTemporalFrequency','9999');
gratings_temporal_frequency = sscanf(str_temp,'%f');
f21_log.gratings_temporal_frequency = gratings_temporal_frequency;

str_temp = f21log_search(lines,'TunningCurve','CRFSize','9999');
crf_size = sscanf(str_temp,'%f');
f21_log.crf_size = crf_size;

str_temp = f21log_search(lines,'TunningCurve','CRFDirection','9999');
crf_dirrection = sscanf(str_temp,'%f');
f21_log.crf_dirrection = crf_dirrection;

str_temp = f21log_search(lines,'TunningCurve','nCRFRatio','9999');
ncrf_ratio= sscanf(str_temp,'%f');
f21_log.ncrf_ratio = ncrf_ratio;

str_temp = f21log_search(lines,'TunningCurve','nCRFSize','9999');
ncrf_size= sscanf(str_temp,'%f');
f21_log.ncrf_size = ncrf_size;

str_temp = f21log_search(lines,'TunningCurve','nCRFDirection','9999');
ncrf_direction= sscanf(str_temp,'%f');
f21_log.ncrf_direction = ncrf_direction;

str_temp = f21log_search(lines,'TunningCurve','nCRFContrast','9999');
ncrf_contrast = sscanf(str_temp,'%f');
log.ncrf_contrast = ncrf_contrast;

str_temp = f21log_search(lines,'TunningCurve','OuterDirection','9999');
outer_direction = sscanf(str_temp,'%f');
f21_log.outer_direction = outer_direction;

str_temp = f21log_search(lines,'TunningCurve','SingleTestTime','1');
single_test_time = sscanf(str_temp,'%d');
f21_log.single_test_time = single_test_time;

str_temp = f21log_search(lines,'TunningCurve','Interval','0');
interval = sscanf(str_temp,'%d');
f21_log.interval = interval;

str_temp = f21log_search(lines,'TunningCurve','DataPoints','1');
data_points = sscanf(str_temp,'%d');
f21_log.data_points = data_points;

str_temp = f21log_search(lines,'TunningCurve','Repeats','1');
repeats = sscanf(str_temp,'%d');
f21_log.repeats = repeats;

str_temp = f21log_search(lines,'TunningCurve','DirectionList','1 0');
list_temp = sscanf(str_temp,'%f');
direction_list = list_temp(2:length(list_temp));
f21_log.direction_list = direction_list;

str_temp = f21log_search(lines,'TunningCurve','SpatialList','1 0');
list_temp = sscanf(str_temp,'%f');
spatial_list = list_temp(2:length(list_temp));
f21_log.spatial_list = spatial_list;

str_temp = f21log_search(lines,'TunningCurve','TemporalList','1 0');
list_temp = sscanf(str_temp,'%f');
temporal_list = list_temp(2:length(list_temp));
f21_log.temporal_list = temporal_list;

str_temp = f21log_search(lines,'TunningCurve','TestSequence','1 0');
list_temp = sscanf(str_temp,'%f');
test_sequence = list_temp(2:length(list_temp));
f21_log.test_sequence = reshape(test_sequence,f21_log.data_points, f21_log.repeats);



sort_sequence = zeros(data_points, repeats);
for i=1:repeats
    tmp = test_sequence(((i-1)*data_points+1):(i*data_points));
    [seq inx] = sort(tmp);
    sort_sequence(:,i) = inx;
end
f21_log.sort_sequence = sort_sequence;

str_temp = f21log_search(lines,'TunningCurve','TestValueSequence','1 0');
list_temp = sscanf(str_temp,'%f');
test_value_sequence = list_temp(2:length(list_temp));
f21_log.test_value_sequence = test_value_sequence;

sort_value_sequence = test_value_sequence(sort_sequence(:,1));
f21_log.sort_value_sequence = sort_value_sequence;

return;



%%

%read f21 log file information into a cell.
%     lines = flog_load(filename);
%
%     INPUTS
%     filename   - filename (string) of the log file
%  
%     OUTPUTS
%     lines - line array 
%
%     if fails, return 0

function lines = f21log_get(filename)

% open file
[fid,message] = fopen(filename,'rt');
if fid < 0
   fprintf('\n%s\n',message);
   lines = 0;
   return;
end

% load file into str
str = fscanf(fid,'%c');

% close file
fclose(fid);

% convert str to line cell array
lines = str2cell(str);

return



function value = f21log_search(lines,section,keyword,default,prompt_not_found);
% FSEARCH_LOG               Retrieve information from f21 log lines
% 
%     value = f21log_search(lines,section,keyword,default,prompt_not_found);
%
%     INPUTS
%     lines      - log line array (return value of fload_log)
%     section    - section name (string) in log file (like [XXX])
%     keyword    - entry name (string) in log file
%     defualt    - return value (string) if expected section-keyword pair not found
%  
%     OUTPUTS
%     value - string returned 
%
%     Yuxi 07.18.2000
%
%     $ Version 1.0 - Yuxi 07.18.2000 - initial version $

if nargin < 4
   fprintf('\nFSEARCH_LOG, no enough input arguments, type ''help fsearch_log'' for help\n\n');
   value = default;
   return;
end

if nargin < 5
   prompt_unfound = 0;
else
   prompt_unfound = prompt_not_found;
end

% length of keyword segment in log file
div = 30;

% flag and counter
count = 1;

% convert section name 'XXX' to '[XXX]' 
section_str = sprintf('[%s]',section);

% not found message
not_found_msg = sprintf('\n  %s\n',['fsearch_log, ' section ', ' keyword ', not found !']);

% get length line cell array
numlines = length(lines);

% Loop through the cell array until section_str and keyword are found
while 1
   % check if exceed lines dimension
   if count > numlines
      if prompt_unfound fprintf(not_found_msg); end
      value = sprintf('%s',default);
      return;
   end
   
   % search section
   if isempty(lines{count})
      % Do nothing if blank line
      count = count + 1;
   elseif strcmp(deblank(lines{count}),section_str)
      count = count + 1;
      while 1
         if count > numlines
            if prompt_unfound fprintf(not_found_msg); end
            value = sprintf('%s',default);
            return;
         end
         
         % search keyword
         if isempty(lines{count})
            count = count +1;
         elseif strcmp(lines{count}(1),'[')
            break;
         elseif strcmp(deblank(lines{count}(1:div)),keyword)
            value = lines{count}(div+1:length(lines{count}));
       %     value = strjust(value,'left');
            value = deblank(value);
            return;
         else
            count = count + 1;
         end
      end
   else
      count = count + 1; 
   end
end

return
% STR2CELL               Convert a string to a cell array of lines.
%
%     C = STR2CELL( STR ) creates a cell array C where each cell contains
%     a line of the string STR.
%
%     C = STR2CELL( STR, OPTS ), where OPTS is 'L', 'T' or both, removes
%     leading and/or trailing blank lines from the string before
%     converting to a cell array.
%
%     If the string contains LFs (linefeed characters, decimal 10), the
%     input string is split at their position after all CRs (carriage
%     return characters, decimal 13) have been removed. If there are no
%     LFs, the string is split at the position of the CRs. This should
%     ensure that the string is split correctly with both UNIX (LF), DOS
%     (CR+LF) and MAC (CR) definitions of a newline.

%     Author:      Peter J. Acklam
%     Time-stamp:  1998-06-30 21:20:01
%     E-mail:      jacklam@math.uio.no
%     WWW URL:     http://www.math.uio.no/~jacklam

function c = str2cell( varargin )

error( nargchk( 1, 3, nargin ) );

%
% Assign default values to parameters that can be changed by command
% line options.
%
strip_lead  = 0;
strip_trail = 0;

%
% Process command line options.
%
while length( varargin ) > 1
   opt = varargin{2};
   if ~ischar( opt )
      error( 'Options must be strings.' );
   end
   switch opt
      case { 'l', 'L' }
         strip_lead = 1;
      case { 'u', 'U' }
         strip_trail = 1;
      otherwise
         error( [ 'Unknown option: ' opt ] );
   end
   varargin(2) = [];
end
str = varargin{1};

%
% Strip leading blank lines.
%
if strip_lead
   k = find( ~isspace(str) );
   if ~isempty(k)
      k = min(k);
      str = str(k:end);
   end
end

%
% Strip trailing blank lines.
%
if strip_trail
   k = find( ~isspace(str) );
   if ~isempty(k)
      k = max(k);
      str = str(1:k);
   end
end

%
% Quick exit if string is empty.
%
if isempty( str )
   c = { '' };
   return
end

%
% Find the characters that separate the lines.
%
k = find( str == 10 );                  % Find all LF (/n) chars. /n
l = find( str == 13 );                  % Find all CR (return) chars.
if isempty( k )                         % If no LF chars were found
   k = l;                               %   split at CR chars.
else                                    % Or else
   if ~isempty( l )                     % If there are CR chars
      str(l) = [];                      % remove them.
      k = find( str == 10 );            % refind all LF chars.
   end
end

%
% Avoid empty last string in output list when string ends in a newline.
%
if ~isempty( k ) & k(end) == length(str)
   k = [ 0 k ];                         % Add beginning.
else
   k = [ 0 k length(str)+1 ];           % Add beginning and end.
end

%
% Now split the string into lines.
%
n = length(k) - 1;                      % Number of lines.
c = cell( n, 1 );                       % Initialize output.
for i = 1:n
   c{i} = str( k(i)+1 : k(i+1)-1 );     % Extract line.
end


