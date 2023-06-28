function [data3,data3Std,angleLeft2]=avgRepeatNO_updated(data2,angleLeft,input)

% This function is called by GaussianFit, which is called by step 4c
    % This fxn gives you average and standard deviation of dff0 values for
        % each stimulus angle, averaged over stim repetitions, for each
        % frame
    % data2 is dff0 values sorted by the frame indices of the diff stimuli
        % in their original sequence
    % angleLeft2 is a vector with sorted stim angles
[m,n,p]=size(data2); %m is total # frames in session (after deleting what you don't want); n=1; p=1
if p>1 % dealing with images -- TYPICALLY p=1, SO THIS IS NOT RELEVANT
    c=unique(angleLeft); %c lists the stim angles w/o repeat
    p2=length(c); %p2 = # stim angles

    data3=zeros(m,n,p2,'uint16'); %will have avg values
    data3Std=zeros(m,n,p2); %will have st error
    for ii=1:p2 %for each stim angle
        cI=(angleLeft==c(ii)); %find indices for this stim angle
        img=nanmean(double(data2(:,:,cI)),3); %data2 is sorted dff0 values
        img2=nanstd(double(data2(:,:,cI)),0,3);
        data3(:,:,ii)=uint16(img);
        data3Std(:,:,ii)=img2;
        if nargin>2
            if input.SEflag==1
                data3Std(:,:,ii)=data3Std(:,:,ii)/sqrt(sum(cI));% I used length(cI) before, which was severely wrong. (R's note)
            end
        end        
    end 
    angleLeft2=c;    
else % dealing with curves (e.g. dff0 values)
    c=unique(angleLeft); %returns the stim angles, w/o repeats
    p2=length(c); %number of stimulus angles used
    data3=zeros(p2,n); %p2 is # of frames in session after deleting unwanted frames; n=1
    data3Std=zeros(p2,n);
    for ii=1:p2  %cycle through all stim angles
        cI=(angleLeft==c(ii)); %set a given angle; find indices in angleLeft corresponding to that angle
        img=nanmean(double(data2(cI,:)),1);  %data2 is dff0 values sorted by frame indices of stimuli
                                             %so you're averaging the dff0
                                             %values for each frame, for a
                                             %given stim angle; you avg
                                             %across frames and across stim
                                             %repeats to get a single value
        img2=nanstd(double(data2(cI,:)),0,1);   %standard error
        data3(ii,:)=img;
        data3Std(ii,:)=img2;
        if nargin>2
            if input.SEflag==1
                data3Std(ii,:)=data3Std(ii,:)/sqrt(sum(cI));
            end
        end         
    end
    angleLeft2=c; 
    
end

