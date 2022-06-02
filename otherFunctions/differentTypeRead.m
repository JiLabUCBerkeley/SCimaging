function img=differentTypeRead(file,fileType)
if nargin==1
    fileType='others';
end
if strcmp(fileType,'txt')
    img=load(file);
elseif strcmp(fileType,'mat')
    img2=load(file);
    img2Name=fieldnames(img2);
    img=img2.(img2Name{1});
else
    img=(imread(file));
%     img=double(imread(file));
%     if size(img,3)~=1
%         img=img(:,:,1);
%     end
end  