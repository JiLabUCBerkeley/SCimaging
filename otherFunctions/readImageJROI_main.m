function [bw,xy]=readImageJROI_main(filePath,file2,img)
[m,n]=size(img);
figure(31);clf;imshow(img,[]);hold on;
fileName=dir(fullfile(filePath,file2));
% deal with single ROI
if isempty(fileName)
    [~,~,ext]=fileparts(file2);
    file2=strrep(file2,ext,'.roi');
    fileName=dir(fullfile(filePath,file2));
    file2=fileName(1).name;
    roi1=ReadImageJROI(fullfile(filePath,file2));
    roi1={roi1};
else
    file2=fileName(1).name;
    roi1=ReadImageJROI(fullfile(filePath,file2));
end

p1=length(roi1);
% p2=2;
xy=cell(p1,1);
% figure(1);clf;imshow(avg,[]);hold on;
% 
bw=false(m,n,p1);
for ii=1:p1
    aaTmp=roi1{ii};
    if strcmp(roi1{ii}.strType,'Polygon') || strcmp(roi1{ii}.strType,'Freehand')
    xy{ii}=roi1{ii}.mnCoordinates;
    
    
    elseif strcmp(roi1{ii}.strType,'Oval')
        yx1=roi1{ii}.vnRectBounds;
        y1=yx1(1);y2=yx1(3);
        x1=yx1(2);x2=yx1(4);
        x0=x1+x2;x0=x0/2;
        y0=y1+y2;y0=y0/2;
        a=x2-x1;a=a/2;
        b=y2-y1;b=b/2;
        theta=linspace(0,2*pi,100);
        theta2=theta(:);
        xy{ii}=[a*cos(theta2)+x0,b*sin(theta2)+y0];
    elseif strcmp(roi1{ii}.strType,'Rectangle') && ~isfield(aaTmp,'vfShapes')
        yx1=roi1{ii}.vnRectBounds;
        y1=yx1(1);y2=yx1(3);
        x1=yx1(2);x2=yx1(4);       
        xx=[x1,x1,x2,x2,x1].';yy=[y1,y2,y2,y1,y1].';
        xy{ii}=[xx,yy];
    elseif strcmp(roi1{ii}.strType,'Rectangle') && isfield(aaTmp,'vfShapes')
        % note: it is not the perfect. we are still missing a small
        % section. 
        a1=aaTmp.vfShapes;
        %% only delete the first 4
        a=[a1;99];        
%         a(a==4)=[];
        ap=length(a);
        b=reshape(a,3,ap/3);
        find4Flag=find(b(1,:)==4);
        b(1,find4Flag(1))=0.001;
        a=b(:);
        a(a==0.001)=[];
        a(end)=[];
        a(end)=[];
        
        ap=length(a);
        b=reshape(a,3,ap/3);        
        xy{ii}=b(2:3,:).';
        
    end
    bw(:,:,ii)=roipoly(img,round(xy{ii}(:,1)),round(xy{ii}(:,2)));
%     bw(:,:,ii)=roipoly(img,xy{ii}(:,1),xy{ii}(:,2));
    plot(xy{ii}(:,1),xy{ii}(:,2));
end
%% get rid of emplty roi
for ii=p1:-1:1
    bw_tmp=bw(:,:,ii);
    if sum(bw_tmp(:))<1
        bw(:,:,ii)=[];
        xy(ii)=[];
    end
end
