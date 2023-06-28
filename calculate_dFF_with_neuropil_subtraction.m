function fileSave=calculate_dFF_with_neuropil_subtraction(filePath,ROI_file,tif_file,options)

% This function calculates dF/F0 for ROIs, using neuropil subtraction.
% 


if nargin == 0  
    filePath='Y:\Katharine\2022_12_17\m952\combined_s05_s08_s09_s10\stack\s10_ortho'; %filePath containing tif stack
    ROI_file='*.zip'; %name of file where ROIs are stored; e.g. zip file from ImageJ 
    tif_file='20*.tif'; %specify the tif file name here. 

    % Parameters for neuropil subtraction    
    NeuropilSize=30; % um over which neuropil signal will be calculated; about 3 times the diameter of cells
    neuropil_coe=0.7; % coefficient for neuropil subtraction; 0.7 works well
end


% Load tif stack
tif_fileName=dir(fullfile(filePath,tif_file)); %tif file
data=read_file(fullfile(filePath,tif_fileName.name));
nFrames = size(data,3);

% Get info from stimulus .mat file
% qt is frames per stimulus, # repetitions per stim, and # of stim angles
[order,qt]=StimulationSequence_moving_bars(filePath,nFrames);
options.qt=qt;
avg=mean(data,3); %average the intensity values of each pixel across frames
stdImg=std(single(data),0,3); %take standard deviation of intensity values of each pixel across frames

%%
% Calculate dFF0

% Read the ROI file
ROI_fileName=dir(fullfile(filePath,ROI_file)); %file from imageJ with ROIs
% Get the ROI masks (bw) and coordinates (xy)
[bw,xy]=readImageJROI_main(filePath,ROI_fileName.name,avg);
nROI=length(xy);
[XpixNum,YpixNum,nFrames]=size(data);
 
% Pre-allocate some arrays
dff0 = zeros(nFrames,nROI); 
f0s=cell(nROI,1);
baseline=zeros(1,nROI);
Intensity_raw=zeros(nFrames,nROI);

% Get raw intensity for each ROI across frames
disp('Calculating intensity for each ROI...')
for ii=1:nROI   
    ROIpos = xy{ii}; %coordinates for ROI
    % Crop each frame around the ROI position; crop_frames has intensity
    % values for the portion of the frame around the given ROI
    crop_frames = data(max(1,floor(min(ROIpos(:,2)))):min(XpixNum,ceil(max(ROIpos(:,2)))),max(1,floor(min(ROIpos(:,1)))):min(YpixNum,ceil(max(ROIpos(:,1)))),:);
    % single_frame_mask is a logical array with the "footprint" of the ROI
    single_frame_mask = bw(max(1,floor(min(ROIpos(:,2)))):min(XpixNum,ceil(max(ROIpos(:,2)))),max(1,floor(min(ROIpos(:,1)))):min(YpixNum,ceil(max(ROIpos(:,1)))),ii);
    % Calculate the mean intensity of the ROI for each frame
    Intensity_raw(:,ii) = squeeze(mean(squeeze(mean(bsxfun(@times, double(crop_frames), double(single_frame_mask))))))*length(single_frame_mask(:))/sum(single_frame_mask(:));
end
 


%% Neuropil subtraction


% Pre-allocate some arrays
baseline_neuropil = zeros(1,nROI);
Intensity_neuropil_subtr=zeros(nFrames, nROI);
Intensity_neuropil = zeros(nFrames, nROI);
neuropil_baseline_vector = zeros(nFrames, nROI);

% Create neuropil mask, which is the inverse of the combination of all ROI
% masks
bw_all=sum(bw,3); % put all roi masks on a single frame
bw_all(bw_all>0)=1;
neuropil_mask=logical(ones(size(bw_all))-bw_all);

for ii = 1:nROI
    crrtROIxy = xy{ii}; %coordinates of ROI
    % Specify an area around the ROI over which neuropil signal will be
    % calculated
    width = 2*NeuropilSize + max(crrtROIxy(:,1)) - min(crrtROIxy(:,1)) ;
    height = 2*NeuropilSize + max(crrtROIxy(:,2)) - min(crrtROIxy(:,2)) ;
    bkgrdROI = getBackground(avg,bw(:,:,ii),xy(ii),3,width,height);
    bwROI=squeeze(bkgrdROI.bwROIAll);
    bwROI=and(bwROI, neuropil_mask); %combine with neuropil mask

    I = bkgrdROI.I;
    I1=I{1};
    [Iy,Ix]=ind2sub([size(data,1),size(data,2)],I1);
    Iy2=min(Iy):max(Iy);
    Ix2=min(Ix):max(Ix);
    crop_back=data(Iy2,Ix2,:);
    crop_mask=bwROI(Iy2,Ix2);
    Intensity_neuropil=squeeze(mean(squeeze(mean(bsxfun(@times, double(crop_back), double(crop_mask))))))*length(crop_mask(:))/sum(crop_mask(:)); % fluorescence intensity of neuropil
    Intensity_thisROI = Intensity_raw(:,ii);
    
    % Get neuropil baseline
    thisIntensity = Intensity_neuropil';
    [N, X] = hist(thisIntensity, 50); %N is the # of counts in each bin; X is the center of each bin **Don't change to histogram**
    f0_mode = find(N==max(N)); %find the index of the bin with the highest peak (i.e. the most common intensity value)
    if length(f0_mode) > 1, f0_mode = f0_mode(1); end  %if there's more than 1 peak in the distribution, just take the first peak
    half_width = (X(end)-X(1))/(length(X)-1)/2; %half the width of each bin
    f0=find(thisIntensity>X(f0_mode)-half_width&thisIntensity<X(f0_mode)+half_width); 
            %find the indices of "thisIntensity" for which the intensity
            %values are within this range of most common values
    baseline_vector = mean(thisIntensity(f0));  %average the intensity values in this range to get the overall mode

    dff0_neuropil = Intensity_neuropil' - baseline_vector;
    neuropil_baseline_vector(:,ii) = baseline_vector;
    baseline_neuropil(ii) = mean(baseline_vector);

    % subtract deltaF of neuropil from ROI:
    Intensity_neuropil_subtr(:,ii) = Intensity_thisROI - neuropil_coe* dff0_neuropil;

end




%% Use the neuropil-subtracted intensity values to calculate dF/F0

dff0 = zeros(nFrames,nROI);
for ii=1:nROI
    display(['calculate df/f0 for ROI ',num2str(ii)])
    Intensity_thisROI = Intensity_neuropil_subtr(:,ii);
    [N, X] = hist(Intensity_thisROI, 50); %N is the # of counts in each bin; X is the center of each bin **Don't change to histogram**
    f0_mode = find(N==max(N)); %find the index of the bin with the highest peak (i.e. the most common intensity value)
    if length(f0_mode) > 1, f0_mode = f0_mode(1); end  %if there's more than 1 peak in the distribution, just take the first peak
    half_width = (X(end)-X(1))/(length(X)-1)/2; %half the width of each bin
    f0=find(Intensity_thisROI>X(f0_mode)-half_width&Intensity_thisROI<X(f0_mode)+half_width); 
            %find the indices of "thisIntensity" for which the intensity
            %values are within this range of most common values
    baseline = mean(Intensity_thisROI(f0));  %average the intensity values in this range to get the overall mode
    dff0(:,ii) = (Intensity_thisROI - baseline)/baseline*100;  %subtract each intensity value by the mode, divide by mode
end




% Save the data
fileSave='ROI_dff0.mat';

save(fullfile(filePath,fileSave),'Intensity_neuropil_subtr', 'baseline', 'f0s', 'avg', 'stdImg', ...
'dff0','bw','xy','nROI','qt','order','tif_file','Intensity_raw', 'Intensity_neuropil','baseline_neuropil','options','neuropil_coe');
