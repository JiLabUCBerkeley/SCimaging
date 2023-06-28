function calculate_orientation_selectivity(filePath,file,options)


close all

if nargin==0
    filePath='Y:\Katharine\2022_12_17\m952\combined_s05_s08_s09_s10\stack\s08_clockwise';
    file='ROI*.mat';
    
    % Delete frames of gray period when grating/bars aren't shown
    % For our protocol, this is the 1st 14 frames of each trial
    options.deleteN=[14,0];  
    
    % Specify direction of bars' movement relative to their orientation; 
    % set both to 0 if bars or grating moves orthogonally to their
    % orientation
    options.rotate_clockwise = 1;
    options.rotate_counterclockwise = 0;

    options.resp_rois_gray = 1; %if set to 1, responsive non-tuned ROIs will be shown as gray in the tuning map                                                               
    options.negativeMethod=1; % 1: set negative df/f0 values to zero when calculating OSI,DSI,GOSI,GDSI; % 2: shift the curve; others: do nothing
    filterMethod.name=0;% 0: do nothing; 1: gaussian; 2: median; 3: wiener filter; 4: max filter;5: butterworth; 
    options.filterMethod=filterMethod;
    
    % Set upper and lower bounds of parameters for gaussian fitting       
    options.LB=[1 -30 5 0.00001 5 -0.2]; %[amp1 theta sigma1 amp2 sigma2 dc];
    options.UB=[1.5 360 180 1.5 180 .5]; %[amp1 theta sigma1 amp2 sigma2 dc];
    % NOTE: to change sigma1 and sigma2 initial values, go to script
    % curveGaussianFitting_v2 and set the values in line ~51-52
end


% Set criteria for determining if an ROI is orientation selective (OS).
% Our default criteria: max dF/F is at least 10%; passes anova, is fit to
% gaussian
criteria.maxYThresh=10; %threshold to determine weather the cell is responsive; default value=10% (dF/F must reach at least this value for an ROI to be considered responsive
OS.annovaFlag=[0.05,1];% the first value is the threshold for the p-value the 2nd value is the option to choose (1) or not to choose (0) to include p value of annova test;
OS.fitFlag=1; % option to use fitted data (1) or raw df/f (0) to calculate OSI itself; 
OS.goodFitFlag=[0.4, 0.6, 1]; % the first value is the threshold of SSE; the second value is the threshold of RSQ; the last one is option to use (1) or not to use (0) critetia of good fitting to designate the cell is OS;
OS.gOSI_OSI=[0.25, 0, 2];% the first value is the threshold of gOSI; the second value is the threshold of OSI; the last one is to choose gOSI>threshold (1) or OSI>threshold (2)
criteria.OS=OS; %all the OS criteria are gathered in a structure called "criteria"
options.criteria=criteria;

%%
% Create folder with current date to put the results
result = ['ROI_maps_created_' datestr(now,'yyyymmdd') '_fitted'];
resultPath_main=fullfile(filePath,result);
if exist(resultPath_main)==7
else
    mkdir(resultPath_main)
end

% Create subfolders; 'map' will have the OS colormaps; 'params' will have
% files with values for OSI, fitted angle, etc.
result_map='map'; 
resultPath_map=fullfile(resultPath_main,result_map);
if exist(resultPath_map)
else
    mkdir(resultPath_map)
end

result_params='params';  
resultPath_params=fullfile(resultPath_main,result_params);
if exist(resultPath_params)
else
    mkdir(resultPath_params)
end



%%
% Load file from previous step that has neuropil-subtracted dff0 values.
% This file will usually be names "ROI*.mat"
fileName=dir(fullfile(filePath,file)); 
load(fullfile(filePath,fileName.name));



nROI=length(xy); %xy gives coordinates for each ROI
warning('off','stats:statrobustfit:IterationLimit');

% Get info from ROI dff0 file; delete unwanted frames (e.g. gray period at start of each trial)
deleteN=options.deleteN; %frames to be deleted
angleNO=qt(3); %how many stim directions were presented; qt comes from ROI_fj file
trialNO=qt(2); %number of trials for each stim
nFrames=length(order)/angleNO/trialNO; %"order" comes from loading the ROI file; total # frames per stim, before frames are trimmed
order=deleteMN(order,deleteN,nFrames); %get rid of deleted frames
[Y,I]=sort(order(:,1));  %sort the 1st col of "order," which is the values for all stim angles
order(:,2:4)=[(1:length(I)).',Y(:)*1,I(:)]; %replace col 2-4 of "order" with intgers from 1 to total # frames, Y, and I, to account for deleted frames 
framePerSti=length(order)/angleNO/trialNO;  % # frames per stim after unwanted frames are trimmed
qt=[framePerSti,angleNO,trialNO];
options.angleNO=angleNO;
options.trialNO=trialNO;
options.framePerSti=framePerSti;
options.qt=qt;



%%
% Fit each ROI to a double gaussian; calculate OSI, DSI, gOSI, gDSI, etc

% Pre-allocate arrays
fitAng=nan(nROI,1);
pvalue=ones(nROI,1);
tic; 
options.plotN=13;
DSI=zeros(nROI,1);
OSI=zeros(nROI,1);
gDSI=zeros(nROI,1);
gOSI=zeros(nROI,1);
pref_direction=zeros(nROI,1);
FWHM=zeros(nROI,1);
minY=FWHM;
minOptions=minY;
maxY=minY;
dff0_avg=zeros(angleNO,nROI);
RSQ=minY;
SSE=minY;
pvalueMin=0.05;
max_dff0=nan(nROI,1);

for ii=1:nROI  
    disp(['analyzing ROI',num2str(ii,'%03d')])
    dff0_this_ROI=dff0(:,ii); %y2 is dff0 values for every frame (no frames omitted)

    dff0_del_frames=deleteMN(dff0_this_ROI,deleteN,nFrames); %delete frames of static period 
    
    
    % Anova test; input is dff0 values ordered according to stimulus angle; 
    % e.g. first are all the values for 0 deg, then 30 deg, etc.
    % Output is dff0 values averaged across frames and averaged 
    % across stimulus repeats to get one dff0 value  for each stimulus
    pvalue(ii)=pvalueGet(dff0_del_frames(order(:,end)),options); 
%     pvalue(ii) = pvalue_thisROI;

    % Gaussian fitting; create a structure called "output" that
    % has all the calculations for each ROI
    output=GaussianFit_updated(dff0_del_frames,order,options);
    SSE(ii)=output.SSE;
    RSQ(ii)=output.RSQ;
    maxY(ii)=output.maxY;
    minY(ii)=output.minY;
    minOptions(ii)=output.minOptions;%=minOptions;
    DSI_fitted(ii)=output.DSI_ii;
    OSI_fitted(ii)=output.OSI_ii;
    gDSI(ii)=output.gDSI;
    gOSI(ii)=output.gOSI;
    fitAng(ii)=output.theta; %fitted angle from gaussian fit (not vector sum; label is confusing)
    pref_direction(ii)=output.pref_direction/pi*180;  %vector sum from polar plot
    FWHM(ii)=output.FWHM;
    dff0_avg(:,ii)=output.dff0_avg; %dff0 averaged across stim repeats and frames
        % Note that dff0_avg should be same as the output derived from
        % the pvalueGet fxn above; there are slight diff's due to rounding
        % errors
 
    if pvalue(ii)<=pvalueMin %get max_dff0 for ROIs that pass anova
        max_dff0(ii)=max(dff0_avg(:,ii));  %map of max dff0 
    end
    
    % Average dff0 across stimulus repeats and calculate standard error
    options.SEflag=1;
    [dff0_avg_stimRepeats,dff0_std_stimRepeats,~]=avgTrials(dff0_del_frames,order(:,1),framePerSti,angleNO,trialNO,options); %data3 is avg dff0 for each stim

   if ii==1 %do this just for the first ROI
       stimAngles=output.x1; %angles used as stimuli
       interpolatedAngles=output.x2; %vector of integers 0:1:360
       dataFit=zeros(length(interpolatedAngles),nROI+1); %interpolated data
       dataFit(:,1)=interpolatedAngles(:);
       dataUnfit=zeros(length(stimAngles),1+2*nROI);
       dataFit(:,1+ii)=output.z2(:); %result of gaussian fit to data
       dataUnfit(:,2*ii:2*ii+1)=[output.y1(:),output.y1SE(:)]; %non-interpolated data
       dataUnfit(:,1)=stimAngles(:);
       dataAndError=zeros(length(dff0_avg_stimRepeats),1+2*nROI);
       dataAndError(:,1)=(1:length(dff0_avg_stimRepeats)).'; %1st col is indices
       dataAndError(:,2*ii:2*ii+1)=[dff0_avg_stimRepeats(:),dff0_std_stimRepeats(:)]; 
    end

   dataFit(:,1+ii)=output.z2(:);
   dataUnfit(:,2*ii:2*ii+1)=[output.y1(:),output.y1SE(:)];  
   dataAndError(:,2*ii:2*ii+1)=[dff0_avg_stimRepeats(:),dff0_std_stimRepeats(:)];
end

% Save avg dff0 
save(fullfile(filePath, 'dff0_avg.mat'),'dff0_avg');

% For rotated bars protocol, pref_direction needs to be adjusted
if options.rotate_clockwise == 1
    pref_direction = pref_direction + 315;
    pref_direction(pref_direction >= 360) = pref_direction(pref_direction >= 360) - 360;
elseif options.rotate_counterclockwise == 1
    pref_direction = pref_direction + 45;
    pref_direction(pref_direction >= 360) = pref_direction(pref_direction >= 360) - 360;
end



%% Calculate raw OSI and DSI 
angleSeries=((1:(angleNO+1))-1)*(360/angleNO); %0:30:360 (all our stimulus angles)

% Find raw preferred angle; ie stimulus angle giving highest response
prefAngle_raw=zeros(nROI,1);
for ii=1:nROI
    dff0_avg_thisROI=dff0_avg(:,ii); %avg dff0 for each stim angle
    difference=abs(angleSeries-fitAng(ii));
    angle_rawTmp=find(difference==min(difference)); %find the smallest difference btn the fitted angle and one of the angles in the series
    % deal with the condition that angle_rawTmp is 360; set it =0
    for jj=1:length(angle_rawTmp)
        if angle_rawTmp(jj)==angleNO+1
            angle_rawTmp(jj)=1;
        end
    end
    % deal with the condition that angle_rawTmp is exactly in the middle of
    % two angles; choose the larger angle
    if length(angle_rawTmp)==2
        if dff0_avg_thisROI(angle_rawTmp(1))>dff0_avg_thisROI(angle_rawTmp(2))
            angle_rawTmp=angle_rawTmp(1);
        else
            angle_rawTmp=angle_rawTmp(2);
        end
    end
    prefAngle_raw(ii)=angle_rawTmp;
end

% Calculate preferred angle +-pi/2 and +pi
step=2*pi/angleNO;  %step size is 30 deg if you have 12 angles
prefAngTmp=exp(1i*(prefAngle_raw-1)*2*pi/angleNO); 

pAng_plusPi=angle(prefAngTmp*exp(1i*pi));
pAng_plusPi(pAng_plusPi<0)=pAng_plusPi(pAng_plusPi<0)+2*pi;
pAng_plusPi=pAng_plusPi/step+1; %range: 1:angleNO+1(1:13)

pAng_plus_halfPi=angle(prefAngTmp*exp(1i*pi/2));
pAng_plus_halfPi(pAng_plus_halfPi<0)=pAng_plus_halfPi(pAng_plus_halfPi<0)+2*pi;
pAng_plus_halfPi=pAng_plus_halfPi/step+1;        

pAng_minus_halfPi=angle(prefAngTmp*exp(-1i*pi/2));
pAng_minus_halfPi(pAng_minus_halfPi<0)=pAng_minus_halfPi(pAng_minus_halfPi<0)+2*pi;
pAng_minus_halfPi=pAng_minus_halfPi/step+1;           

prefAngle_raw=round(prefAngle_raw);
pAng_plusPi=round(pAng_plusPi);
pAng_plus_halfPi=round(pAng_plus_halfPi);
pAng_minus_halfPi=round(pAng_minus_halfPi);   

% If prefAngle_raw is 360, make it 0 deg
prefAngle_raw(prefAngle_raw==(angleNO+1))=1;
pAng_plusPi(pAng_plusPi==(angleNO+1))=1;
pAng_plus_halfPi(pAng_plus_halfPi==(angleNO+1))=1;
pAng_minus_halfPi(pAng_minus_halfPi==(angleNO+1))=1;  

% Find dff0 at each angle
for ii=1:nROI
    z1=dff0_avg(prefAngle_raw(ii),ii);
    z2=dff0_avg(pAng_plusPi(ii),ii);
    z3=dff0_avg(pAng_plus_halfPi(ii),ii);
    z4=dff0_avg(pAng_minus_halfPi(ii),ii);  
    % Calculate raw OSI
    if z1+z2+z3+z4==0
        OSI_raw(ii)=(z1+z2-z3-z4)/(z1+z2+z3+z4+10*eps);
    else
        OSI_raw(ii)=(z1+z2-z3-z4)/(z1+z2+z3+z4);
    end
    % Calculate raw DSI
    if z1+z2==0
        DSI_raw(ii)=abs(z1-z2)/(z1+z2+10*eps);
    else
        DSI_raw(ii)=abs(z1-z2)/(z1+z2);
    end            
end


%% Determine which ROIs are OS 
% The function OSCriteria_kat creates a structure called OS that contains 2 substructures
%    Substructure "result" is a logical array indicating whether a given
%        ROI is OS (1) or not OS (2)
%    Substructure "output" contains 2 things:
%        (1) Structure "flagHead" will heading labels for (a) isOS? 
%            (b) dF/F>10% (c) pvalue<0.05 (d) good fit? (e) fitted OSI>0?
%        (2) Structure OS_flagAll: this is a logical array with 5 columns,
%             each matching one of the headers and designating whether a
%             given ROI meets the criteria or not
%
handles.maxY=maxY;
handles.dff0_avg=dff0_avg;
handles.pvalue=pvalue;
handles.RSQ=RSQ;
handles.SSE=SSE;
handles.OSI_raw=OSI_raw(:);
handles.OSI_fitted=OSI_fitted(:);
handles.gOSI=gOSI;
handles.shuffle=[];

[OS.result,OS.output]=OSCriteria_moving_bars(handles,options.criteria);

%% Find centroids of ROIs 
centroids = zeros(size(bw,3),2); 
for ii = 1:size(bw,3)
    tmp = regionprops(bw(:,:,ii),'Centroid'); 
    centroids(ii,:) = tmp.Centroid;
end


%% Create Excel spreadsheet and .mat files with all the data
xlsHead={'ROI index', 'fittedAngle', 'pvalue','DSI_fitted','DSI_raw','OSI_fitted','OSI_raw', ...
    'gDSI','gOSI','theta(pref_direction)','FWHM','maxY','minY_avgForFit', ...
    'minOptions(0:N/A;1:shift;2:negSet0)','centroidX','centroidY'};
xlsData=[(1:nROI).',fitAng(:),pvalue(:),DSI_fitted(:),DSI_raw(:),OSI_fitted(:),OSI_raw(:), ...
    gDSI(:),gOSI,pref_direction(:),FWHM(:),maxY(:),minY(:), ...
    minOptions(:),centroids(:,1),centroids(:,2)];
xlsHead=[xlsHead,OS.output.flagHead];
xlsData=[xlsData,OS.output.OS_flagAll];

% Save the data
xlsDataSave=[xlsHead;num2cell(xlsData)];
xlswrite(fullfile(resultPath_params,'allPar.xls'),xlsDataSave); %save Excel file
save(fullfile(resultPath_params,'allPar.mat'),'xlsHead','xlsData','OS'); %save .mat file


%%
% Create colorscale maps of ROIs showing preferred angle for OS and DS
options.resultPath=resultPath_map;
options.fileName=fileName;
pvalueMinFlag=pvalue<=pvalueMin;
OSFlag=OS.result;

%  OSI colorMap
options.fitAng=fitAng;
options.avg=avg;  
options.stdImg=stdImg;
options.pvalueMinFlag=pvalueMinFlag;
options.bw=bw;    
options.OSFlag=OSFlag;
options.respFlag = OS.output.OS_flagAll(:,2);
% options.respFlag=DS.output.DS_flagAll(:,2); %these are all responsive ROIs
colorMapforROIs(options);

  

  
%%
function colorMapforROIs(input)
% This fxn creates various images that show a "map" of the ROIs,
% some of them overlaid on the average or standev image of the
% stack.
fitAng=input.fitAng; %values for fitted angle
avg=input.avg;  %avg intensity of pixels in 2p image
bw=input.bw; %masks for ROIs to show shape

OSFlag=input.OSFlag; %ROIs that are OS
respFlag=input.respFlag; %ROIs that are responsive

% Get indices of OS ROIs
OS_ROIs=find(OSFlag); 
nOS_ROIs=length(OS_ROIs); % # of OS ROIs

% Get indices of nontuned responsive ROIs
if input.resp_rois_gray==1
    resp_ROIs=find(respFlag & ~OSFlag); 
end


% Set up maps, same # pixels as 2p image, that will
% show locations of ROIs
[x_pix,y_pix]=size(avg);
map=nan(x_pix,y_pix);  


% Find positions of OS ROIs
centroidXY=zeros(nOS_ROIs,2);
for ii=1:nOS_ROIs
    jj=OS_ROIs(ii);
    tmp=regionprops(bw(:,:,jj),'centroid');
    centroidXY(ii,:)=tmp.Centroid;
    map(bw(:,:,jj))= fitAng(jj); %ROI masks now contain values for fitted angle
end

% Find positions of gray ROIs (responsive but not OS)
map_gray=map; %start with OS map
if input.resp_rois_gray==1
    for ii=1:length(resp_ROIs)
        jj=resp_ROIs(ii);
        map_gray(bw(:,:,jj))=540;
    end
end



%% Create map of OS ROIs on 0-180 scale 
input.fullColor=1;
input.lhFlag=0;
img=avg*0+255;
img=uint8(img);
[rgbImg,map2]=mapHSV_180Degrees(img,map_gray,input);
fig1 = figure;
imshow(uint8(rgbImg*256));

% Save the OS colormap in various formats
resultPath=input.resultPath;
saveas(fig1,fullfile(resultPath,'OSI_ColorMap.fig'))
exportgraphics(fig1,fullfile(resultPath,'OSI_ColorMap.eps'))
saveas(fig1,fullfile(resultPath,'OSI_ColorMap.pdf'))
saveas(fig1,fullfile(resultPath,'OSI_ColorMap.png'))



