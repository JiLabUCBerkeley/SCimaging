function [map,OSI]=OSI_standard(dff0,point_angle,minOSI)
    %this code is to calculate the tuned angle of cells. OSI is defined as
    %(DF/F_max-DF/F_othorg)/(DF/F_max+DF/F_othorg); DF/F_othorg is defined
    %as angle of max DF/F+-90 degrees. dff0 is mXn matrix. 1st dimension is
    %time and the second dimension is roi or pixel; minOSI is the threshold
    % 
[m,n]=size(dff0); % 1st dimension : time; 2nd dimension: roi or pixel
map=nan(1,n);
point=point_angle(1);
angleNO=point_angle(2);
dff02=reshape(dff0,point,angleNO,n);

[Y]=max(dff02,[],1);
Y=squeeze(Y); %amplitude, 1st dimension: angle; 2nd dimenson: ROIs; 
[maxF,I]=max(Y,[],1);% find maximum DF/F and corresponding angle
% find +-90degrees apart
ii3=I-90/(360/angleNO);

ii3Tmp=ii3<1; % deal with the condition when angle is <0 degree
ii3(ii3Tmp)=ii3(ii3Tmp)+angleNO;

ii2=I+90/(360/angleNO);
ii2Tmp=ii2>angleNO;
ii2(ii2Tmp)=ii2(ii2Tmp)-angleNO;% deal with the condition when angle is > 360 degree
minF=I*0;
% minF2=I*0;
for ii=1:length(I)
    minF1=mean(Y([ii2(ii),ii3(ii)],ii));
    minF(ii)=minF1;
end
maxF=double(maxF);
OSI=double(maxF)-minF;
OSI=OSI./(maxF+minF);
OSI_flag=OSI>=minOSI;
map(OSI_flag)=(I(OSI_flag)-1)*360/angleNO;

