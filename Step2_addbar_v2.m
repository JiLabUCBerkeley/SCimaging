function Step2_addbar_v2(filePath,file,options)
% Step2_addbar_v2 does three jobs: 1) add moving bars to raw tif stack; 2)
% sort images according to stimulus angles, do averaging across trials, and
% then add bars to average images; 3) similar to 2), but add bars to std
% images;
% this code will go through all sessions in the same folder

if nargin==0
    filePath='D:\Ji lab Liang\results\20160710 352411 dLGN 25x\session4_200um\stack'; 
    file='20160710*.tif';
    filePath='C:\Users\lur4\Box\superior colliculus manuscript Liang\data\figure1_tuning';
    file='2018*.tif';    
        options.preMoving=3;% number of static frames before stimuli move
        options.postMoving=0;% number of static frames after stimuli move
        options.blank=0;    % if it is 0, pre and post moving stimuli are static grating; if it is 1, then it is blank;
end
fileName=dir(fullfile(filePath,file)); % get file name
options.noBarSaveFlag=0;
h_wait = waitbar(0.2,'checking image quantities');
p=length(fileName);
imgN=zeros(p,1);
for ii=1:p
    info=imfinfo(fullfile(filePath,fileName(ii).name));
    imgN(ii)=length(info);
end
 waitbar(0.3, h_wait,'get sequence');
%% get sequence
% order is a cell type; length of the order is corresponding to the session
% number;  each element of order has 4 columns. 
% The first column is the angle sequence in degree for individual images; the
% second column is the index (1,2,3,4,... to number of frames); the third column is one by
% sorting the 1st column, and the fourth column is reorganized index of the
% 2nd column; order{1}(:,3)=order{1}(order{1}(:,4),1);
% point_trial_angle is p by 3 matrix; p is the session number; the 3rd column is
% the angle quantity, the 2nd column is the trial repetition number; and
% the 1st column is the number of frames per stimulus.
[order,point_trial_angle]=StimulationSequence(filePath,fileName,imgN);
trial_angle=point_trial_angle(:,2:3);%  get the trial repetition number and the angle quantity
framePerSti=point_trial_angle(:,1); %get the the number of frames per stimulus
%%
p=length(fileName);
for ii=1:p
    waitbar(ii/p, h_wait,[num2str(floor(ii/p*100),'%02d'),' completed']);
    file1=fileName(ii).name;    
    file2=['bar_',file1];
        data=loadImgSequence(filePath,file1);% load data
      if options.noBarSaveFlag==1
          [data_avg,~,angleLeft_tmp]=avgTrials(data,order{ii}(:,1),framePerSti(ii),trial_angle(ii,2),trial_angle(ii,1));
             saveSingleTif(fullfile(filePath,strrep(file2,'bar_','AvgNOBar_')), data_avg)
              [data_avg_avg]=avgRepeatNO(data_avg,angleLeft_tmp);
              saveSingleTif(fullfile(filePath,strrep(file2,'bar_','Avg_AvgNOBar_')), data_avg_avg)
         end       
        % add gratings to raw images
        data=addOrientation_v2(data,order{ii}(:,1),trial_angle(ii,2),framePerSti(ii),options);
        % average across trials, sort images
        % according to their angle of grating
        [data_avg,data_Std,angleLeft]=avgTrials(data,order{ii}(:,1),framePerSti(ii),trial_angle(ii,2),trial_angle(ii,1));
        % add gratings to the averged and sorted images
         data_avg_bar=addOrientation_v2(data_avg,angleLeft,trial_angle(ii,2),framePerSti(ii),options);

         saveSingleTif(fullfile(filePath,strrep(file2,'bar_','AvgBar_')), data_avg_bar)
          data_Std_bar=addOrientation_v2(data_Std,angleLeft,trial_angle(ii,2),framePerSti(ii),options);
          saveSingleTif(fullfile(filePath,strrep(file2,'bar_','StdBar_')), data_Std_bar)
    saveSingleTif(fullfile(filePath,file2), data)
    t=order{ii};
    file2b=['StimulusOrder',file1];
    file2b=strrep(file2b,'.tif','.txt');
    save(fullfile(filePath,file2b),'t','-ASCII');
end
close(h_wait)

