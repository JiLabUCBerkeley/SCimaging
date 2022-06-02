function [rgbImg,map2,discreteMap,discreteMap2]=mapHSV(img,angleMap,options)
ratio=1; % changed on 05/29/2018 to avoid the same color at two ends
stretchFlag_user=0;
    if nargin>2
        if isfield(options,'lhFlag')
            if options.lhFlag==1
                stretchFlag_user=1;
            end
        end
    end
if stretchFlag_user==1
    img2=mat2gray(double(img),[options.low,options.high]);
else
    img2=mat2gray(double(img));
end

%
%% processing discrete map
discreteFlag=0;
    if nargin>2
        if isfield(options,'discreteFlag')
            if options.discreteFlag==1
                discreteFlag=1;
            end
        end
    end
if discreteFlag==0
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
%     tmp=linspace(0,1,n/2);
%     h=ones(m,1)*[tmp,tmp(end:-1:2)];
%     h=ones(m,1)*[tmp,tmp(1:end)];
    h=ones(m,1)*linspace(0,ratio,n);
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
else
    [m,n]=size(img2);
    mapH=options.mapH;
    mapH2=[mapH(:);1];
    angleMap2=angleMap;
    
    % angleMap2=angleMap;
    % [m1,n1]=size(angleMap2);

    angleMap2(isnan(angleMap))=length(mapH2);
    angleMap2=mapH2(angleMap2);
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

        m=10;n=360;

    h=ones(m,1)*linspace(0,ratio,n);
    [m,n]=size(h);
    s=ones(m,n);
    v=ones(m,n)*1;
    map2=zeros(m,n,3);
    map2(:,:,1)=h;
    map2(:,:,2)=s;
    map2(:,:,3)=v;
    map2=hsv2rgb(map2);
    

   k=30;
%      m1=3;n1=3; n1=round(length(mapH(:))/m1);
   for ss=3:5
       if mod(length(mapH(:)),ss)==0
           m1=ss;
           n1=length(mapH(:))/m1;
           break;
       end
   end

    h=reshape(mapH(:).',m1,n1);
    h=imresize(h,k,'nearest');
    
    [m,n]=size(h);
    s=ones(m,n);
    v=ones(m,n)*1;
    discreteMap=zeros(m,n,3);
    discreteMap(:,:,1)=h;
    discreteMap(:,:,2)=s;
    discreteMap(:,:,3)=v;
    discreteMap=hsv2rgb(discreteMap);    
%%
m1=1;n1=length(mapH(:))/m1;
    h=reshape(mapH(:).',m1,n1);
    h=imresize(h,k,'nearest');
    
    [m,n]=size(h);
    s=ones(m,n);
    v=ones(m,n)*1;
    discreteMap2=zeros(m,n,3);
    discreteMap2(:,:,1)=h;
    discreteMap2(:,:,2)=s;
    discreteMap2(:,:,3)=v;
    discreteMap2=hsv2rgb(discreteMap2);   
end