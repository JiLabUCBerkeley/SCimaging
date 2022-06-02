function [resultPath,resultName]=step1_convertToTiffStack_Unregistered(filePath,file)
% load individual tif files, register the stack and then save it
% image stack without registration is also saved
% average image and std image are also saved
% copy any log files to the same folder
if nargin==0
    filePath='C:\serverData\2016.08.01 blood vesselSC';
    file='20160801_16*.tif';
    % set parameters for registration; the ierations stops when either
    % condition is met. 

end
tic;
data=single(loadImgSequence(filePath,file));% load individual images;
t=toc;
disp(['loading images spent ',num2str(t,'%.0f'),'seconds'])
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

resultName=strrep(resultName,'.tif','unReg.tif');

saveSingleTif(fullfile(resultPath,resultName), data);
dataMean=mean(data,3);
% combine unregistered image and registered image side by side;

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
%copy log files
fileName2=dir(fullfile(filePath,'*.log'));
if ~isempty(fileName2)
%     p=length(fileName2);
    copyfile(fullfile(filePath,'*.log'),resultPath);
end