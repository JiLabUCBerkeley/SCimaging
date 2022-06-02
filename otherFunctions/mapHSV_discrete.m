function [rgbImg,map2]=mapHSV_discrete(img,angleMap,input)
    stretchFlag_user=0;
    if nargin>2
        if isfield(input,'lhFlag')
            if input.lhFlag==1
                stretchFlag_user=1;
            end
        end
    end
if stretchFlag_user==1
    img2=mat2gray(double(img),[input.low,input.high]);
else
    img2=mat2gray(double(img));
end

[m,n]=size(img2);

angleMap2=mod(angleMap,180*2)/(180*2);
% angleMap2=angleMap;
% [m1,n1]=size(angleMap2);

angleMap2(isnan(angleMap))=1;
h=angleMap2;
s=ones(m,n);
s(isnan(angleMap))=0;
% s(:,end)=1;
% h(:,end)=linspace(0,1,m);
v=img2;
% if isfield(input,'brighterFlag')
%     t=isnan(angleMap);
%     t2=v(t);
%     t2(t2<0.2)=t2(t2<0.2)*2;
%     v(t)=t2;
% end

hsvImg=zeros(m,n,3);
hsvImg(:,:,1)=h;
hsvImg(:,:,2)=s;
hsvImg(:,:,3)=v;
rgbImg=hsv2rgb(hsvImg);
figure(11);clf;imshow(rgbImg);

m=10;n=362;
tmp=linspace(0,1,n/2);
h=ones(m,1)*[tmp,tmp(end:-1:2)];
h=ones(m,1)*[tmp,tmp(1:end)];
h=ones(m,1)*linspace(0,1,n);
[m,n]=size(h);
s=ones(m,n);
v=ones(m,n)*1;
map2=zeros(m,n,3);
map2(:,:,1)=h;
map2(:,:,2)=s;
map2(:,:,3)=v;
map2=hsv2rgb(map2);

% [X,map]=rgb2ind(rgbImg,128);
% subplot(1,2,2);imshow(X,map);colorbar;
% trash=0;