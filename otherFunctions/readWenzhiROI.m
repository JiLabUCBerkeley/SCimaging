function [bw,xy]=readWenzhiROI(filePath,file,img)
% modifed on August 16 2016, to correct the issue that cannot process
% single point ROI
fileName=dir(fullfile(filePath,[file,'*.mat']));
file1=fileName(1).name;
matObj=matfile(fullfile(filePath,file1));
ROIs=matObj.ROIs;
roi1=ROIs.ROIinfo;
%% for testing
% roi1=roi1(1:20);
[m,n]=size(img);
figure(31);clf;imshow(img,[]);hold on;
p1=length(roi1);
xy=cell(p1,1);
bw=false(m,n,p1);
for ii=1:p1
    xy{ii}=roi1(ii).ROIpos;
    if size(xy{ii},1)==1 && size(xy{ii},2)==2
        xyTmp=round(xy{ii});
        bw(xyTmp(2),xyTmp(1),ii)=true;
        plot(xy{ii}(:,1),xy{ii}(:,2),'+');        
    else
        bw(:,:,ii)=roipoly(img,xy{ii}(:,1),xy{ii}(:,2));
        plot(xy{ii}(:,1),xy{ii}(:,2));        
    end

end
drawnow;