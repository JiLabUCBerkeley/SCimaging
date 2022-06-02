function [resultPath,resultName]=step1_convertToTiffStack_registered(filePath,file,options)
% load individual tif files, register the stack and then save it
% image stack without registration is also saved
% average image and std image are also saved
% copy any log files to the same folder
if nargin==0
    filePath='\\dm11\genie\Yajie_GENIE_stuff\project\backupdatafromoldprojects\20150813\20150813 anm313917 592um ccneurons';
    file='2015*.tif';
%     filePath='C:\serverData\2015.05.11 b'; file='05_Ori_D093*.tif';
    % set parameters for registration; the ierations stops when either
    % condition is met. 
    options.maxLoop=15; % maximum iterations allowed to register images; 
    options.minSSE=5; % threshold of SSE; if SSE is smaller than this value, then ierations stop;    
    options.previousRefFlag=0; % when it is 0, register each image with respect to the previous image; when it is 1, then register each image with respect to average of all images;
end
ta=datevec(now);
data=single(loadImgSequence(filePath,file));% load individual images;
tb=datevec(now);
disp(['loading images spent ',num2str(etime(tb,ta),'%.0f'),'seconds'])
%% create folder name
result='stack';
resultPath=fullfile(filePath,result);
if exist(resultPath)==7
else
    mkdir(resultPath)
end

fileName=dir(fullfile(filePath,file));% get file names;
fileTmp=fileName(1).name;
resultNameTmp=fileTmp;
% truncate the name, for example:
%20160123_01_D90_L20_BCr1c4_1.2x_Visible_000001_VisStim_ada_rep01_01_2016-01-23_181058_000_INT1 _XY
% throw away _000001_VisStim_ada_rep01_01_2016-01-23_181058_000_INT1 _XY
id0=regexp(resultNameTmp,'_\d\d\d\d\d\d_'); 
if isempty(id0)
    resultName=resultNameTmp;
else
    resultName=[fileTmp(1:id0-1),'.tif'];
end


% register images
tic;
% for ii=1:10
%     data(:,:,ii)=circshift(data(:,:,ii),25);
% end

dataRaw=data;
% []
[m,n,p]=size(data);
if p>10000
    [data,yxShiftAll]=imgRegistration_largeDataSet(data,options);
else
    [data,yxShiftAll]=imgRegistration(data,options);
end

% [data,yxShiftAll]=imgRegistration(data,options);
t=toc;
disp(['registering images spent ',num2str(t,'%.0f'),'seconds'])

% save regiserted images
saveSingleTif(fullfile(resultPath,resultName), data);
save(fullfile(resultPath,['shift_',strrep(resultName,'.tif','.txt')]),'yxShiftAll','-ASCII');%20190118 added'-v7.3'
dataMean=mean(data,3);
dataMax=max(dataMean(:));
% combine unregistered image and registered image side by side;

saveCmprFlag=0;%% was 1 20180712
if isfield(options,'noSaveUnregistered')
    if options.noSaveUnregistered==1;
        saveCmprFlag=0;
    end
end
if saveCmprFlag==1
    dataAll=zeros(m,2*n+1,p,'single');
    dataAll(:,1:n,:)=dataRaw;

    dataAll(:,n+1,:)=dataMax;
    dataAll(:,n+2:end,:)=data;
    saveSingleTif(fullfile(resultPath,['cmpr_',resultName]), dataAll);
end

%% save average image of all registered images
result='avg';
resultPathavg=fullfile(resultPath,result);
if exist(resultPathavg)==7
else
    mkdir(resultPathavg)
end
saveSingleTif(fullfile(resultPathavg,resultName),dataMean);

%% save standard deviation image of all registered images
result='std';
resultPathstd=fullfile(resultPath,result);
if exist(resultPathstd)==7
else
    mkdir(resultPathstd)
end
tic
saveSingleTif(fullfile(resultPathstd,resultName), std(data,0,3));

t=toc;
disp(['saving images spent ',num2str(t,'%d'),'seconds'])

% try
%copy log files
fileName2=dir(fullfile(filePath,'*.log'));
if ~isempty(fileName2)
%     p=length(fileName2);
    copyfile(fullfile(filePath,'*.log'),resultPath);

end
% catch me
% end
% try

%copy MAT files
fileName2=dir(fullfile(filePath,'*.mat'));
if ~isempty(fileName2)
%     p=length(fileName2);
    copyfile(fullfile(filePath,'*.mat'),resultPath);

end
% catch me
% end








