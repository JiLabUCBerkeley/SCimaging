function [reg_frames,yxShiftAll]=imgRegistration(frames2reg,options)

% To calcium imaging, how to select a good template is critical for the motion
% correction by using cross correlation. Here we use mean of whole stack as
% template to loop the calculation until all motion is corrected.

% Created by Wenzhi Sun, Aug 26 2013.
% Changed to a fast calculation by Uri Dubin, Feb 13 2014.
% Modified to removed shifted edges by Wenzhi Sun Apr 10 2014.
% Remove the crop code, 04/05/2015. To remove the artifact of circ_shift,
% the images need to be embbeded before registration.

% dtfreg_fast_mean - like original algrism but with FFT tricks and fast max find

% Input:
% frames2reg - images to be corrected
% template - outer source template
% Output:
% yxShift - Nx2 shift array in x and y
% reg_frames - registered image frames. same size as input stack.

%  if nargin < 2, iftrim = 0; end

% modifed by Rongwen Lu on 07/29/2016
% 1 fix the bug that can only handle images with square size;
% 2 add options for the maxLoop; default value is 10;

% modifed by Rongwen Lu on 08/11/2016
% 1 add option to register currrent image agains previous one, useful for
% 3D structrual data

% modified by Raphael Turcotte on 04/06/2017
% add function imgRegistration_extRef()
% Perform the registration using an external reference

if nargin==1
    options.maxLoop=10;
    options.previousRefFlag=0;
    options.extRef=0;
elseif nargin==2
    if ~isfield(options,'previousRefFlag')
        options.previousRefFlag=0;
    end
    if ~isfield(options,'extRef') || ~isfield(options,'extRefImage')
        options.extRef=0;
    end
end

if size(frames2reg,3)==1
    reg_frames=frames2reg;
    yxShiftAll=[0 0];
    return; 
end

if isfield(options,'smallFOV')
else
    options.smallFOV.flag=0;
end

smallFOV=options.smallFOV;

if options.smallFOV.flag==1
   if options.previousRefFlag==1
        upsampling=1;
        [reg_frames,yxShiftAll]=imageRegistration_previousRef(frames2reg(smallFOV.y,smallFOV.x,:),upsampling);
    elseif options.previousRefFlag==0
        [reg_frames,yxShiftAll]=imgRegistration_avgRef(frames2reg(smallFOV.y,smallFOV.x,:),options);
    end
else
    if options.previousRefFlag==1
        upsampling=1;
        [reg_frames,yxShiftAll]=imageRegistration_previousRef(frames2reg,upsampling);
    elseif options.previousRefFlag==0
        if options.extRef == 0
            [reg_frames,yxShiftAll]=imgRegistration_avgRef(frames2reg,options);
        elseif options.extRef == 1
            [reg_frames,yxShiftAll]=imgRegistration_extRef(frames2reg,options);
        end
    end
    return;
end



if isfield(options,'smallFOV')
    [m,n,p]=size(frames2reg);
    reg_frames=zeros(m,n,p,'single');
    if options.smallFOV.flag==1
        
        for m=1:size(reg_frames,3)
            reg_frames(:,:,m)=circshift(frames2reg(:,:,m),yxShiftAll(m,:));
        end
        
    end
end

function [reg_frames,yxShiftAll]=imgRegistration_avgRef(frames2reg,options)

if ~isfield(options,'minSSE')
    options.minSSE=0;
end
[~,~,nT]=size(frames2reg);
yxShift=zeros(nT,2);
yxShiftAll=zeros(nT,2);
reg_frames=single(frames2reg);

% transform entire data array; To speed up, this tranformatin is done only
% once.
frames2reg_fft=fft2(reg_frames);

SSE = 10^4;

t=1;
while SSE>options.minSSE && t<=options.maxLoop
    %template for registration
    template_fft=conj(fft2(nanmean(reg_frames,3)));
    
    this_yxShift=dtfreg_fast_shifts(frames2reg_fft,template_fft);
    for m=1:nT
        reg_frames(:,:,m)=circshift(frames2reg(:,:,m),this_yxShift(m,:));
    end
    
    SSE=sum(sum((this_yxShift-yxShift).*(this_yxShift-yxShift)));
    disp(['iteration  ',num2str(t),';SSE=',num2str(SSE)]);
    yxShift = this_yxShift;
    yxShiftAll=this_yxShift;
    t=t+1;
end

reg_frames = cast(reg_frames,'like',frames2reg);

function [yxShift]= dtfreg_fast_shifts(frames2reg_fft,template_fft,varargin)

% dtfreg_mean_Fast: to calculate the y-x shifts for each frames to the template with DTF
% Input:
% frames2reg_fft - fft transformed images to be registered
% template - FFTed template
% Output:
% yxShift - N x 2 shift array in x and y
% modified on 01/02/2015

[nC,nR,nT]=size(frames2reg_fft);
[m,n,o]=size(frames2reg_fft);
% yxShift=zeros(nT,2);

% multiply by template
frames2reg_fft=bsxfun(@times,frames2reg_fft,template_fft);
frames2reg=frames2reg_fft;
% tic;
% parfor ii=1:nT
% %     frames2reg(:,:,ii)=
%     frames2reg(:,:,ii)=ifft2(frames2reg_fft(:,:,ii),'symmetric');
% end
% toc
% tic
frames2reg=ifft2(frames2reg_fft,'symmetric');
% toc
if nargin>=3
    rongwenCorrected=1;
else
    rongwenCorrected=0;
end
if nC==nR && rongwenCorrected==0
    [max1,loc1]=max(frames2reg,[],1);
    max1=squeeze(max1);
    loc1=squeeze(loc1);
    
    [max2,loc2]=max(max1,[],1);
    rloc=loc1(sub2ind([nC nT],loc2,1:nT));
    cloc=loc2;
    
    rloc=rloc-1;
    cloc=cloc-1;
    iBool=rloc>fix(nR/2);
    rloc(iBool)=rloc(iBool)-nR; %   row_shift = rloc - nR - 1;
    iBool=cloc>fix(nC/2);
    cloc(iBool)=cloc(iBool)-nC; %   col_shift = cloc - nC - 1;
    
    yxShift=-round([rloc(:) cloc(:)]);
else
    yxShift=zeros(o,2);
    %% to fix Wenzhi's code's issue that only can register images with the nC=nR
    %     CC = ifft2(buf1ft.*conj(buf2ft));
    for ii=1:o
        CC=frames2reg(:,:,ii);
        [max1,loc1] = max(CC);
        [max2,loc2] = max(max1);
        rloc=loc1(loc2);
        cloc=loc2;
        md2 = fix(m/2);
        nd2 = fix(n/2);
        if rloc > md2
            row_shift = rloc - m - 1;
        else
            row_shift = rloc - 1;
        end
        
        if cloc > nd2
            col_shift = cloc - n - 1;
        else
            col_shift = cloc - 1;
        end
        yxShift(ii,:)=-[row_shift,col_shift];
        
    end
    
end
% find maxima (3d dim is time)

function [reg_frames,yxShiftAll]=imgRegistration_extRef(frames2reg,options)

if ~isfield(options,'minSSE')
    options.minSSE=0;
end
[~,~,nT]=size(frames2reg);
yxShift=zeros(nT,2);
yxShiftAll=zeros(nT,2);
reg_frames=single(frames2reg);

% transform entire data array; To speed up, this tranformatin is done only
% once.
frames2reg_fft=fft2(reg_frames);

SSE = 10^4;

t=1;
while SSE>options.minSSE && t<=options.maxLoop
    %template for registration
    template_fft=conj(fft2(options.extRefImage));
    
    this_yxShift=dtfreg_fast_shifts(frames2reg_fft,template_fft);
    for m=1:nT
        reg_frames(:,:,m)=circshift(frames2reg(:,:,m),this_yxShift(m,:));
    end
    
    SSE=sum(sum((this_yxShift-yxShift).*(this_yxShift-yxShift)));
    disp(['iteration  ',num2str(t),';SSE=',num2str(SSE)]);
    yxShift = this_yxShift;
    yxShiftAll=this_yxShift;
    t=t+1;
end

reg_frames = cast(reg_frames,'like',frames2reg);
