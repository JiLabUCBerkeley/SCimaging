% Combine data from various sessions into one file
% You first need to copy every allPar file that you want to include into the filePath

filePath = ''; % Directory where data files are stored
files = dir([filePath '\allPar*.mat']); %this creates a structure with all the files in the filePath
nFiles = size(files,1); %how many allPar files are in the folder
data = [];

% Load each file and append all the data to the xlsData_combROIs array
for i = 1:nFiles
    load(fullfile(filePath,files(i).name), 'xlsData');
    data = [data; xlsData];
end
load(fullfile(filePath,files(1).name), 'xlsHead'); %get the column labels
save(fullfile(filePath, 'allPar_vert.mat'),'data','xlsHead'); %save as a new file


