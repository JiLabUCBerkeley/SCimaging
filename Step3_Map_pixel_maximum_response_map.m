function Step3Map_pixel_maximum_response_map(filePath0,file0,options)
% generate max response map, set the scale bar at line 17
if nargin==0
        filePath0='G:\superior colliculus manuscript Liang\figures\Figure 1\fig1d_rawfiles\rawdata_generatingOSmap'; 
    file0='2018*.tif';
    filterMethod.name=0;% 0: do nothing; 1: gaussian; 2: median; 3: wiener
    % filter; 4: max fitler;5: butterworth; 
    filterMethod.name=0;filterMethod.sigmaNumber=1;filterMethod.sizeNumber=3;
    options.filterMethod=filterMethod;
    options.deleteN=[3,0];% remove first 2 and last 0
    options.nanROI=nan; % set nan 
     bck.method=2; % 1: mode; 2: fixed, e.g., gray period, it uses the first number of options.deleteN to decide the gray period. 3: 10% mimimum;4:to be determined;
     bck.grayNumber=options.deleteN(1);
     options.bck=bck;
%      grayN=100;%%line 145 change the number for gray color map

     options.maximumHist=[0:5:50];%%line 145 change the number for gray color map
     options.refreshFlag=1;
     
     options.pvalueMin=1;
     options.minDff=10;
end
filePath{1}=filePath0;
file{1}=file0;

p=length(filePath);
for ii=1:p
    if ~isempty(filePath{ii})
        fileName=dir(fullfile(filePath{ii},file{ii}));
        p1=length(fileName);
        for jj=1:p1
            orientationMap_generalized_sub(filePath{ii},fileName(jj),options);
        end
    end
end

function orientationMap_generalized_sub(filePath,fileName,options)
    id=1:11;
    pvalueMin=options.pvalueMin;
if options.refreshFlag==1
    data=loadImgSequence(filePath,fileName(1).name);
    [m,n,p]=size(data);
    [order,point_trial_angle]=StimulationSequence(filePath,fileName,p);
    filterMethod=options.filterMethod;
    data=differentTypeReadFilter(data,filterMethod);
    order=order{1};
    options.kk=1;%% spines
    avg=nanmean(data,3);
    stdImg=nanstd(single(data),0,3);
   
    %  [a, MSGID] = lastwarn();
    warning('off','stats:statrobustfit:IterationLimit');
    bw1=true(m,n);
    % method 1, mod; most frequent mod
    %method 2 fixed: fixPersti has two inputs; fixPersti(1) is the static frame number per stimulus;
    % fixPersti(2) is the moving frame number per stimulus;
    % method 3: fix is a vector, defing the frame numbers as baseline; example:
    % fix=[1,10,11,12], then average the frame # 1, 10, 11, 12 as baseline
    % method 4:  use the minimum 20% values as baseline;
    method=1;%mod
    deleteN=options.deleteN;
    nanROI=options.nanROI;
    angleNO=point_trial_angle(3);
    trialNO=point_trial_angle(2);
    framePerSti2=length(order)/angleNO/trialNO;
    order=deleteMN(order,deleteN,framePerSti2);
    [Y,I]=sort(order(:,1));
    order(:,2:4)=[(1:length(I)).',Y(:)*1,I(:)]; 
    framePerSti=length(order)/angleNO/trialNO;
    point_trial_angle=[framePerSti,angleNO,trialNO];
    options.angleNO=angleNO;
    options.trialNO=trialNO;
    options.framePerSti=framePerSti;
    options.qt=point_trial_angle;
    
    map=nan(m,n);
    pvalue=ones(m,n);
    tic;
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        c=parcluster;
        parpool('local',c.NumWorkers);
    else
    %     poolsize = poolobj.NumWorkers
    end
    % for ii=1:m
    bck=options.bck;
    % pvalue=zeros();

    parfor ii=1:m
        bck2=bck;
        disp(['analyzing line',num2str(ii,'%03d')])
        for jj=1:n
            if bw1(ii,jj)
                % obtain temporal curve for individual pixels;
                y=squeeze(data(ii,jj,:));% reduce the redundant dimension
                y2=double(y(:));
    %             y2=calculatedDff0(y2,method);% change to dff0
            if bck2.method==1
                y2=calculatedDff0(y2,method);% change to dff0
            elseif bck2.method==2
                bck2.fixPersti=[bck2.grayNumber,framePerSti2];
                y2=calculatedDff0(y2,bck2.method,bck2);% change to dff0
            elseif bck2.method==3
                y2=calculatedDff0(y2,4);%
            elseif bck2.method==5

            end
    %             brob=robustfit(imgM,y2);

    %             y2=y2-(brob(1)+brob(2)*imgM);

                if length(nanROI)>1
                    y2(nanROI,:)=nan;
                else
                end

                y3=deleteMN(y2,deleteN,framePerSti2);
    %             y3=;
    % get pvalue for anova test
                [pvalue(ii,jj),y4]=pvalueGet(y3(order(:,end)),options); 
                % do double Gaussian fitting if pvalue is small than pvalueMin
%                 if pvalue(ii,jj)<=pvalueMin
    %                 y4=reshape(y,framePerSti,trialNO,angleNO);
    %                 output=GaussianFit(y3,order,input);
                    map(ii,jj)=max(y4);
%                 end
            end


        end
    end
     fileSave=['maxMap_',strrep(fileName(1).name,'.tif','.mat')];
    save(fullfile(filePath,fileSave),'map','avg','stdImg','filterMethod','pvalue');   
else
   fileName2=dir(fullfile(filePath,['maxMap_',fileName(1).name(id),'*.mat']));
    load(fullfile(filePath,fileName2(1).name));  
end
% [m,n]=size(map);
map(pvalue>pvalueMin)=nan;
map(map<=options.minDff)=nan;
map2_save=map;
figure;%%20180907 added the following four lines
imagesc(map2_save);
colormap(gray);
colorbar;
caxis([0 50]);
cmbFile1=['gray_map2',fileName(1).name(id),'.fig'];
saveas(gcf, fullfile(filePath,cmbFile1));

% file1=strrep(file,'map.tif','pValue_color.mat');
% save(fullfile(resultPath,['max_',file1]),'map','pvalue','rgbImg','stdImg2');
options.mapH=([43,233,177,311,58,118,1,280,80]-1)/360;options.discreteFlag=1;
maximumHist=options.maximumHist;

p=length(maximumHist);
for ii=1:p-1
    if ii==1
        map(map<=maximumHist(1))=ii;
    end
    map(map>=maximumHist(ii) & map<maximumHist(ii+1))=ii;
    if ii==p-1
        map(map>=maximumHist(p))=ii;
    end
end
p=p-1;
if length(options.mapH)<p
    options.mapH=linspace(0,300,p)/360;
end
%% processing avg
[m,n]=size(avg);
avg=mat2gray(avg);
avgTmp=sort(avg(:));
ratioAvg=0.99;
avgUpperLimit=avgTmp(round(ratioAvg*m*n));
avg=mat2gray(avg,[0,avgUpperLimit]);
[rgbImg,map2,discreteMap,discreteMap2]=mapHSV(avg,map,options);
[stdImg2,map2,discreteMap,discreteMap2]=mapHSV(stdImg,map,options);
% stdImg2=mapHSV(stdImg,map);
% filterMethod=options.filterMethod;
if filterMethod.name==0
    result='pixelMap';
else
    result=['pixelMap_flted',num2str(filterMethod.name)];
end
if options.deleteN(1)>0
    result=[result,'_trimp',num2str(options.deleteN(1))];
end
resultPath=fullfile(filePath,result);
if exist(resultPath)
else
    mkdir(resultPath)
end
id=1:11;
%%
file=[fileName(1).name(id),'_colormap.tif'];
cmbFile=['max_cmb_',fileName(1).name(id),'.tif'];
imwrite(uint8(rgbImg*255),fullfile(resultPath,cmbFile),'tif','compression','none');
imwrite(uint8(map2_save),fullfile(resultPath,'gray_map2.tif'),'tif','compression','none');



imwrite(uint8(discreteMap2*255),fullfile(resultPath,'max_map2.tif'),'tif','compression','none');
imwrite(uint8(discreteMap*255),fullfile(resultPath,'max_map1.tif'),'tif','compression','none');
% file=fullfile(resultPath,file);
imwrite(uint8(stdImg2*255),fullfile(resultPath,['max_',strrep(file,'.tif','_std.tif')]),'tif','compression','none');



%%
img=avg*0+255;
img=uint8(img);

[rgbImg]=mapHSV(img,map,options);


imwrite(uint8(rgbImg*255),fullfile(resultPath,['max_',strrep(file,'.tif','_360Fullmap.tif')]),'tif','compression','none');
% map2=map(:);
map2_save(isnan(map2_save))=[];
figure(5);clf;histogram(map2_save);
figure(6)
c=gray;
colormap(c)

trash=0;



