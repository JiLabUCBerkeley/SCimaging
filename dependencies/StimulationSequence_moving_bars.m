function [order,point_trial_angle,output]=StimulationSequence_moving_bars(filePath,nFrames)

% This function obtains the order in which visual stimuli are presented.

% "order" is a cell array; length of the order is corresponding to the session
% number;  each element of order has 4 columns. 
% The first column is the sequence of stimulus angles in the order they were presented; the
% second column is the index (1,2,3,4,... to number of frames); 3rd col is
% the sorted stimulus angles, 4th col is the indices corresponding to
% sorted values

stim_file=dir(fullfile(filePath,'*.mat')); %mat file generated by psychophysics toolbox
logFiles=fullfile(filePath,stim_file(1).name);

% Get info from stimulus file
matObj=matfile(logFiles);
output=load(logFiles);

orderTmp1=matObj.angle_sequence; %get stimuli sequence from stumulus file
orderTmp2=orderTmp1(:); 
trial_angle=[matObj.nRepeats,matObj.nAngles]; %trial is # repeats for each stim, angle is # of unique stim
framesPerTrial=nFrames./trial_angle(:,1)./trial_angle(:,2); 
framesPerTrial=floor(framesPerTrial);
orderTmp3=ones(framesPerTrial,1)*orderTmp2.'; %each trial is now a column with #rows = framesPerTrial
orderTmp4=orderTmp3(:);

[Y,I]=sort(orderTmp4); %Y is the angle values; I the indices
order=[orderTmp4(:)*1,(1:length(I)).',Y(:)*1,I(:)]; 
point_trial_angle=[framesPerTrial,trial_angle];




