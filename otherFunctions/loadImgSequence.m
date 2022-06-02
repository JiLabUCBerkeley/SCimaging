function [imgAll]=loadImgSequence(filePath,file)


fileName=dir(fullfile(filePath,file));

p=length(fileName);
infoAll=cell(p,1);
p2=zeros(p,1);
for ii=1:p
    file=fileName(ii).name;
    infoAll{ii}=imfinfo(fullfile(filePath,file));
    p2(ii)=numel(infoAll{ii});
    if ii==1
        if p2(1)==1
    p2=ones(p,1);
    break;

        end
        
        
        
        
        
    end
end
pp=sum(p2);
m=infoAll{1}(1).Height;n=infoAll{1}(1).Width;
avg=zeros(m,n);


imgAll=zeros(m,n,pp,'uint16');

% [bw,xy]=readImageJROI_main(filePath2,file2,avg);

% nROI=length(xy);
% data=zeros(pp,nROI);
% bw_i=cell(nROI,1);
% for ii=1:nROI
%     bw1=bw(:,:,ii);
%     bw_i{ii}=find(bw1(:)==true);
% end

p3=[0;p2];
for ii=1:p
    tic;
    img=readSingleTif(fullfile(filePath,fileName(ii).name));
    y1=sum(p3(1:ii))+1;y2=y1:y1+p3(ii+1)-1;
    imgAll(:,:,y2)=uint16(img);
    t=toc;
    if mod(ii,200)==0
        disp(['loading file',num2str(ii),'using time',num2str(t)])
    end
    
    %     trash=0;
    
    
    
    
end


