function Step4_c_ROIanalysis_99prct_visualdriven_fixOSInegative(filePath0,file0,options)
%% modifed on 08/12/2016, 
% 1) add an option to have Wenzhi's ROIs
% 2) create tuning map by putting color of tuning orientations for individual ROIs
% 3) add ROI number on top of tuning map
% 4) add tuning angles on top of tuning map
% modifed on 08/15/2015,
% correct for the batch processing error;
%modified on 2019023
%line99: increase to 6 frames

if nargin==0
    filePath0='C:\Users\yajie.liang\Box\superior colliculus manuscript Liang\figures\Figure 1\fig1d_rawfiles\tuningmapandcurve';
    file0='ROI_fj*.mat';
    filterMethod.name=0;% 0: do nothing; 1: gaussian; 2: median; 3: wiener
% filter; 4: max fitler;5: butterworth; 
    options.filterMethod=filterMethod;
    options.deleteN=[3,0];% remove first 2 and last 0
    options.nanROI=nan; % set nan 
    options.maximumHist=[0:5:250];
end

options.epsSave=1;
options.curveSave=1;
%% negative value issue
% error: not used at all
options.negativeMethod=1; % 1: set negative value as zero when calculating OSI,DSI,GOSI,GDSI; % 2: shift the curve; others: do nothing
%% boundary        
options.LB=[1 -30 5 0.00001 5 -0.2]; %[amp1 theta sigma1 amp2 sigma2 dc];
options.UB=[1.5 360 180 1.5 180 .5]; %[amp1 theta sigma1 amp2 sigma2 dc];
%         LB = [1 -30 5 0.00001 5 0]; % origonal value
%         UB = [1.5 360 180 1.5 180 .2];origonal value

options.DSI_thresh=0.5;
    fiji_filePath=filePath0; %ROIs drawn in imagej should be saved here with the name"20160707_03*"
%% OS
criteria.maxYThresh=10; % threshold to determine weather the cell is responsive; default value: 10 (10%)
OS.annovaFlag=[0.05,1];% the first value is the threshold for the p-value the 2nd value is the option to choose (1) or not to choose (0) to include p value of annova test;
OS.fitFlag=1; % option to use fitted data (1) or raw df/f (0) to caluclate OSI itself; 
OS.goodFitFlag=[0.4, 0.6, 0]; % the first value is the threshold of SSE; the second value is the threshold of RSQ; the last one is option to use (1) or not to use (0) critetia of good fitting to designate the cell is OS;
OS.gOSI_OSI=[0.25, 0, 2];% the first value is the threshold of gOSI; the second value is the threshold of OSI; the last one is to choose gOSI>0.25 (1) or OSI>0.5 (2)
criteria.OS=OS;

%% DS

DS.annovaFlag=[0.05,1];% the first value is the threshold for the p-value the 2nd value is the option to choDSe (1) or not to choDSe (0) to include p value of annova test;
DS.fitFlag=1; % option to use fitted data (1) or raw df/f (0) to caluclate DSI itself; 
DS.goodFitFlag=[0.4, 0.6, 0]; % the first value is the threshold of SSE; the second value is the threshold of RSQ; the last one is option to use (1) or not to use (0) critetia of good fitting to designate the cell is DS;
DS.gDSI_DSI=[0.25, 0.5, 2];% the first value is the threshold of gDSI; the second value is the threshold of DSI; the last one is to choDSe gDSI>0.25 (1) or DSI>0.5 (2)
DS.OSFlag=0; % if it is 0, ignore whether the cell is OS or not; if it is 1, then only count those who are OS;
criteria.DS=DS;
options.criteria=criteria;
if isfield(options,'DSI_thresh')
else
    options.DSI_thresh=0.5;
end
filePath{1}=filePath0;
file{1}=file0;
options.filePath2{1}=fiji_filePath;
p=length(filePath);
for ii=1:p
    if ~isempty(filePath{ii})
%         matToTifstack_barAdded(filePath{ii},strrep(file{ii},'*_noBar.tif','*.mat'));
        fileName=dir(fullfile(filePath{ii},file{ii}));
        p1=length(fileName);
        for jj=1:p1
            options.jj=jj;
            options.filePath2{jj}=filePath{ii};
            orientationMap_generalized_sub(filePath{ii},fileName(jj),options);
        end
    end
end
end

function orientationMap_generalized_sub(filePath,fileName,options)
disp([filePath,'start----------------------------------------------------------'])
% id=1:11;
fileName2=dir(fullfile(filePath,fileName(1).name));
load(fullfile(filePath,fileName2(1).name));
file=fileName(1).name;
id=1:6;
file(id)=[];
file=strrep(file,'.mat','.tif');
id=1:11;
fileName(1).name=file;
%% for testing
% xy(21:end)=[];

p1=length(xy);
warning('off','stats:statrobustfit:IterationLimit');
%mod
deleteN=options.deleteN;
nanROI=options.nanROI;
angleNO=qt(3);
trialNO=qt(2);
framePerSti2=length(order)/angleNO/trialNO;

order_save=order;

order_deletemoving=movingdeleteMN(order_save,deleteN,framePerSti2);%%20171204
order_deletestation=stationdeleteMN(order_save,deleteN,framePerSti2);%%20171204 the first 3 frames in the moving session were used;%20190123, the first 6 frames

order=deleteMN(order,deleteN,framePerSti2);
[Y,I]=sort(order(:,1));
order(:,2:4)=[(1:length(I)).',Y(:)*1,I(:)]; 


[Y,I]=sort(order_deletemoving(:,1));%%frames only with stationary frames
order_deletemoving(:,2:4)=[(1:length(I)).',Y(:)*1,I(:)]; 
[Y,I]=sort(order_deletestation(:,1));%%frames only with the first 3 frames in the moving session
order_deletestation(:,2:4)=[(1:length(I)).',Y(:)*1,I(:)]; 


framePerSti=length(order)/angleNO/trialNO;
framePerSti_deletemoving=length(order_deletemoving)/angleNO/trialNO;
framePerSti_deletestation=length(order_deletestation)/angleNO/trialNO;

 evokepvalue_paired=zeros(p1,12);  
 evokepvalue_nonpaired=zeros(p1,12);

qt=[framePerSti,angleNO,trialNO];
options.angleNO=angleNO;
options.trialNO=trialNO;
options.framePerSti=framePerSti;
options.qt=qt;
% pvalueMin=0.01;
map=nan(p1,1);
pvalue=ones(p1,1);
visualdrivenpvalue=ones(p1,1); %%20171202
tic;
nanRatio=1;
options.plotN=13;
DSI=zeros(p1,1);
OSI=zeros(p1,1);
gDSI=zeros(p1,1);
gOSI=zeros(p1,1);
pref_direction=zeros(p1,1);
FWHM=zeros(p1,1);
minY=FWHM;
minOptions=minY;
maxY=minY;
dff0_avg=zeros(angleNO,p1);
RSQ=minY;
SSE=minY;
pvalueMin=0.05;
map_max=nan(p1,1);
morethan99pct=zeros(p1,1);%%20171202
evokedpoint5=zeros(p1,1);
evokedpoint1=zeros(p1,1);
evoked_nonpairedpoint5=zeros(p1,1);
evoked_nonpairedpoint1=zeros(p1,1);

for ii=1:p1
    disp(['analyzing ROI',num2str(ii,'%03d')])
    if ii==224
        trash=0;
    end
    y2=dff0(:,ii);
%     y=Intensity(:,ii);
    
            if length(nanROI)>1
                y2(nanROI,:)=nan;
            else
            end
             y3=deleteMN(y2,deleteN,framePerSti2);
             [pvalue(ii),y4]=pvalueGet(y3(order(:,end)),options); 
             
          %% to get the visual evoked p value   
              y3_stationdelete=stationdeleteMN(y2,deleteN,framePerSti2);
              y3_movingdelete=movingdeleteMN(y2,deleteN,framePerSti2);
    
              AA=y3_stationdelete(order_deletestation(:,end));

              y44=reshape(AA,framePerSti_deletestation,trialNO,angleNO);
              y5=squeeze(mean(y44,1));% first dimension is trialNO and the 2nd dimension is the angleNO
    
              BB=y3_movingdelete(order_deletemoving(:,end));
    
              y6=reshape(BB,framePerSti_deletemoving,trialNO,angleNO);
              y7=squeeze(mean(y6,1));% first dimension is trialNO and the 2nd dimension is the angleNO
    
              pp1=angleNO;
    
    
                for iii=1:pp1
                [h,p]=ttest(y5(:,iii),y7(:,iii));
                evokepvalue_paired(ii,iii)=p;
                [h,p]=ttest2(y5(:,iii),y7(:,iii));
                evokepvalue_nonpaired(ii,iii)=p;
                
                end    
            
             aaaa=evokepvalue_paired(ii,:)<0.05;
             evokedpoint5(ii)=sum(aaaa,2);
             bbbb=evokepvalue_paired(ii,:)<0.01;
             evokedpoint1(ii)=sum(bbbb,2);
             
             aaaa=evokepvalue_nonpaired(ii,:)<0.05;
             evoked_nonpairedpoint5(ii)=sum(aaaa,2);
             bbbb=evokepvalue_nonpaired(ii,:)<0.01;
             evoked_nonpairedpoint1(ii)=sum(bbbb,2);
                                            
             
             
             %%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
             output=GaussianFit(y3,order,options);
             SSE(ii)=output.SSE;
             RSQ(ii)=output.RSQ;
             maxY(ii)=output.maxY;
             minY(ii)=output.minY;%=min(ydata);
            minOptions(ii)=output.minOptions;%=minOptions;
             DSI(ii)=output.DSI_ii;
             OSI(ii)=output.OSI_ii;
             gDSI(ii)=output.gDSI;
             gOSI(ii)=output.gOSI;
             map(ii)=output.theta;
             pref_direction(ii)=output.pref_direction/pi*180;
             FWHM(ii)=output.FWHM;
             dff0_avg(:,ii)=output.dff0_avg;
                if pvalue(ii)<=pvalueMin

                    map_max(ii)=max(y4);
                end
                             
                dff0raw=dff0(:,ii);%%%20171202
                 aaa=mean(dff0raw(find(dff0raw>prctile(dff0raw,99))));
                 morethan99pct(ii)=aaa;
    


             %%
              options.SEflag=1;
              [data3,data3SE,angleLeft2]=avgTrials(y3,order(:,1),framePerSti,angleNO,trialNO,options);
              data3b=insertNan(data3,framePerSti,nanRatio);
              data3SEb=insertNan(data3SE,framePerSti,nanRatio);
              %%

             if ii==1
                 x1=output.x1;
                 x2=output.x2;
                 dataFit=zeros(length(x2),p1+1);
                 dataFit(:,1)=x2(:);
                 dataUnfit=zeros(length(x1),1+2*p1);
                 dataFit(:,1+ii)=output.z2(:);
                 dataUnfit(:,2*ii:2*ii+1)=[output.y1(:),output.y1SE(:)];
                 dataUnfit(:,1)=x1(:);
                 
                 dataAndError=zeros(length(data3b),1+2*p1);
                 dataAndError(:,1)=(1:length(data3b)).';
                 dataAndError(:,2*ii:2*ii+1)=[data3b(:),data3SEb(:)];
             end
                 dataFit(:,1+ii)=output.z2(:);
                 dataUnfit(:,2*ii:2*ii+1)=[output.y1(:),output.y1SE(:)];  
                 dataAndError(:,2*ii:2*ii+1)=[data3b(:),data3SEb(:)];
end
[m,n]=size(avg);
map_max2=nan(m,n);    
for ii=1:p1
    if isnan(map_max(ii))
    else
%         tmp=regionprops(bw(:,:,ii),'centroid');
        map_max2(bw(:,:,ii))= map_max(ii);
    end
end
map_max=map_max2;
map_max_save=map_max;
map_max_save(isnan(map_max_save))=[];
figure(34);clf;histogram(map_max_save(:));
options.mapH=([43,233,177,311,58,118,1,280,80]-1)/360;options.discreteFlag=1;
maximumHist=options.maximumHist;
p=length(maximumHist);
for ii=1:p-1
    if ii==1
        map_max(map_max<=maximumHist(1))=ii;
    end
    map_max(map_max>=maximumHist(ii) & map_max<maximumHist(ii+1))=ii;
    if ii==p-1
        map_max(map_max>=maximumHist(p))=ii;
    end
end
p=p-1;
if length(options.mapH)<p
    options.mapH=linspace(0,300,p)/360;
end


OSI_fitted=OSI;
DSI_fitted=DSI;
%% analysis on OS cells
handles.maxY=maxY;
handles.dff0_avg=dff0_avg;
handles.pvalue=pvalue;
handles.RSQ=RSQ;
handles.SSE=SSE;

    x=((1:(angleNO+1))-1)*(360/angleNO);
    preferredAngle_raw=zeros(p1,1);
    for ii=1:p1
        y_avg=dff0_avg(:,ii);
        xTrash=abs(x-map(ii));
%         xTrash2=[xTrash;xTrash(1)];
%% note: map is the theta angle, which is based on the fitting. Sometimes fitting is not very good, which makes the theta inaccurate. 
        angle_rawTmp=find(xTrash==min(xTrash));
        
        angle_rawTmp=find(y_avg==max(y_avg));
        
        %% deal with the condition that angle_rawTmp is 360;
        for jj=1:length(angle_rawTmp)
            if angle_rawTmp(jj)==angleNO+1
                angle_rawTmp(jj)=1;
            end
        end
        %% deal with the condition that angle_rawTmp is exactly in the middle of two points
        if length(angle_rawTmp)==2
            
            if y_avg(angle_rawTmp(1))>y_avg(angle_rawTmp(2))
                angle_rawTmp=angle_rawTmp(1);
            else
                angle_rawTmp=angle_rawTmp(2);
            end
           
        end
         preferredAngle_raw(ii)=angle_rawTmp;

        
    end
step=2*pi/angleNO;
y1=preferredAngle_raw;
y1Tmp=exp(1i*(y1-1)*2*pi/angleNO);
y2=angle(y1Tmp*exp(1i*pi));
y2(y2<0)=y2(y2<0)+2*pi;
y2=y2/step+1;

y3=angle(y1Tmp*exp(1i*pi/2));
y3(y3<0)=y3(y3<0)+2*pi;
y3=y3/step+1;        

y4=angle(y1Tmp*exp(-1i*pi/2));
y4(y4<0)=y4(y4<0)+2*pi;
y4=y4/step+1;           

y1=round(y1);y2=round(y2);y3=round(y3);y4=round(y4);
trash=0;    

y1(y1==(angleNO+1))=1;
y2(y2==(angleNO+1))=1;
y3(y3==(angleNO+1))=1;
y3(y3==(angleNO+1))=1;

for ii=1:p1
%             disp(num2str(ii))
    z1=dff0_avg(y1(ii),ii);

    z2=dff0_avg(y2(ii),ii);
    z3=dff0_avg(y3(ii),ii);
    z4=dff0_avg(y4(ii),ii);   
    if z1+z2+z3+z4==0
        OSI(ii)=(z1+z2-z3-z4)/(z1+z2+z3+z4+10*eps);
    else
        OSI(ii)=(z1+z2-z3-z4)/(z1+z2+z3+z4);
    end
    if z1+z2==0
        DSI(ii)=abs(z1-z2)/(z1+z2+10*eps);
    else
        DSI(ii)=abs(z1-z2)/(z1+z2);
    end            
%             DSI=(z(1)-z(2))/(z(1)+z(2));

end
%% OS
OSI_raw=OSI;
handles.OSI_raw=OSI_raw;
handles.OSI_fitted=OSI_fitted;
handles.gOSI=gOSI;
% criteria=options.criteria;
% OS=criteria.OS;
[OS.result,OS.output]=OSCriteria(handles,options.criteria);

%% DS
DSI_raw=DSI;
handles.DSI_raw=DSI_raw;
handles.DSI_fitted=DSI_fitted;
handles.gDSI=gDSI;
handles.OS=OS;
% criteria=options.criteria;
% DS=criteria.DS;
[DS.result,DS.output]=DSCriteria(handles,options.criteria);

%%
result='curve';
% result='fiji_curves_12plot_overlay';
resultPath1=fullfile(filePath,result);
if exist(resultPath1)==7
else
    mkdir(resultPath1)

end
result='map';
resultPath=fullfile(resultPath1,result);
if exist(resultPath)
else
    mkdir(resultPath)
end

result='params';
resultPath2=fullfile(resultPath1,result);
if exist(resultPath2)
else
    mkdir(resultPath2)
end
%%
avg2=avg;
[m,n]=size(avg2);
avg2=mat2gray(avg2);
avgTmp=sort(avg2(:));
ratioAvg=0.99;
avgUpperLimit=avgTmp(round(ratioAvg*m*n));
avg2=mat2gray(avg2,[0,avgUpperLimit]);
[rgbImg,map2,discreteMap,discreteMap2]=mapHSV(avg2,map_max,options);
[stdImg2]=mapHSV(stdImg,map_max,options);
imwrite(uint8(rgbImg*255),fullfile(resultPath,['max_',fileName(1).name(id),'avg.tif']),'tif','compression','none');
imwrite(uint8(stdImg2*255),fullfile(resultPath,['max_',fileName(1).name(id),'std.tif']),'tif','compression','none');
imwrite(uint8(discreteMap*255),fullfile(resultPath,['max_map1.tif']),'tif','compression','none');
imwrite(uint8(discreteMap2*255),fullfile(resultPath,['max_map2.tif']),'tif','compression','none');
img=avg2*0+255;
img=uint8(img);

[rgbImg]=mapHSV(img,map_max,options);
imwrite(uint8(rgbImg*255),fullfile(resultPath,['max_',fileName(1).name(id),'white.tif']),'tif','compression','none');

options.discreteFlag=0;
options.resultPath=resultPath2;
options.fileName=fileName;
options.pvalue=pvalue;
options.DSI=DSI;
options.OSI=OSI;
options.gDSI=gDSI;

fig=51;

options.plot_condition=1;

%% DSI
options.condition=DS.output.DS_flagAll(:,2); % responsive cells.
options.conditionSplit=DS.result; % DS cells
if options.criteria.DS.fitFlag==1
    options.dataSave=DSI_fitted;
else
    options.dataSave=DSI_raw;
end

options.prefix='DSI';options.fig=fig;fig=fig+1;
output_tmp=drawParameters(options);
DS.dataSplit=output_tmp.dataSplit;
DS.dataRest=output_tmp.dataRest;
DS.dataResponsive=output_tmp.dataResponsive;

options.dataSave=gDSI;options.prefix='gDSI';options.fig=fig;fig=fig+1;
gDSI_save=drawParameters(options);



%% OSI
options.condition=DS.output.DS_flagAll(:,2); % responsive
options.conditionSplit=OS.result;

if options.criteria.OS.fitFlag==1
    options.dataSave=OSI_fitted;
else
    options.dataSave=OSI_raw;
end

options.prefix='OSI';options.fig=fig;fig=fig+1;
output_tmp=drawParameters(options);
OS.dataSplit=output_tmp.dataSplit;
OS.dataRest=output_tmp.dataRest;
OS.dataResponsive=output_tmp.dataResponsive;
options.dataSave=gOSI;options.prefix='gOSI';options.fig=fig;fig=fig+1;
gOSI_save=drawParameters(options);



options.dataSave=pref_direction;options.prefix='theta';options.fig=fig;fig=fig+1;
drawParameters(options);



options.dataSave=pvalue;options.prefix='pvalue';options.fig=fig;fig=fig+1;
drawParameters(options);

%% FWHM
options.condition=DS.output.DS_flagAll(:,2); % responsive
options.dataSave=FWHM;options.prefix='FWHM';options.fig=fig;fig=fig+1;
drawParameters(options);

%% FWHM split into OS and non-OS subsets
options.condition=DS.output.DS_flagAll(:,2); % responsive
options.conditionSplit=OS.result;
options.dataSave=FWHM;options.prefix='FWHMsub';options.fig=fig;fig=fig+1;
FWHM_save=drawParameters(options);
%% fited angle
options.condition=DS.output.DS_flagAll(:,2); % responsive
options.dataSave=map;options.prefix='fitAng';options.fig=fig;fig=fig+1;
drawParameters(options);
%% fited angle, divided into DS, OS excluding DS
options.condition=DS.result; % responsive
options.dataSave=map;options.prefix='fitAngDS';options.fig=fig;fig=fig+1;
drawParameters(options);

OS_excluding_DS=OS.result;
OS_excluding_DS(DS.result)=false;
options.condition=OS_excluding_DS; % responsive
options.dataSave=map;options.prefix='fitAngOS_noDS';options.fig=fig;fig=fig+1;
drawParameters(options);

%% fitted angle, AS and DS cells
AS_and_DS=OS_excluding_DS|DS.result;
options.condition=AS_and_DS;
options.dataSave=map;options.prefix='fitAng_AS_and_DS';options.fig=fig;fig=fig+1;
drawParameters(options);


options.dataSave=maxY(:);options.prefix='maxDF_F';options.fig=fig;fig=fig+1;
drawParameters(options);

[mm1,nn1] =size(evokepvalue_paired);
headSub=cell(1,nn1);
for ii=1:nn1
    headSub{ii}=['pvalue',num2str(ii)];
end
xlsHead={'ROI index', 'fittedAngle', 'pvalue','DSI_fitted','DSI_raw','OSI_fitted','OSI_raw', ...
    'gDSI','gOSI','theta(pref_direction)','FWHM','maxY','minY_avgForFit', 'mean99prctile','evokepoint5','evokedpoint1', 'evoked_nonpairedpoint5','evoked_nonpairedpoint1',...
    'minOptions(0:N/A;1:shift;2:negSet0)'};
xlsHead=[xlsHead,headSub];
xlsData=[(1:p1).',map(:),pvalue(:),DSI_fitted(:),DSI_raw(:),OSI_fitted(:),OSI_raw(:), ...
    gDSI(:),gOSI,pref_direction(:),FWHM(:),maxY(:),minY(:), morethan99pct(:),evokedpoint5(:), evokedpoint1(:),evoked_nonpairedpoint5(:),evoked_nonpairedpoint1(:), ...
    minOptions(:),evokepvalue_paired];
xlsHead=[xlsHead,OS.output.flagHead];
xlsData=[xlsData,OS.output.OS_flagAll];
xlsHead=[xlsHead,DS.output.flagHead];
xlsData=[xlsData,DS.output.DS_flagAll];

% output.flagHead=flagHead;
% output.OS_flagAll=OS_flagAll;

xlsDataSave=[xlsHead;num2cell(xlsData)];

id=1:11;
file=fileName(1).name(id);
xlswrite(fullfile(resultPath2,['allPar_',file,'.xls']),xlsDataSave);
save(fullfile(resultPath2,['allPar_',file,'.mat']),'xlsHead','xlsData','DS','OS','gOSI_save','gDSI_save','FWHM_save');
% filePathTmp=cd;
% cd(resultPath2)
% save ['allPar_',file] 
% cd(filePathTmp)
options.resultPath=resultPath;
% drawOSI(input);
% drawgDSI(input);
%%



pvalueMinFlag=pvalue<=pvalueMin;
OSFlag=OS.result;

options.map=map;
options.avg=avg;  
options.stdImg=stdImg;
options.pvalueMinFlag=pvalueMinFlag;
options.bw=bw;    
options.OSFlag=OSFlag;



colorMapforROIs(options);

%%  DSI colorMap
DSI_flag=DSI<=options.DSI_thresh;
DSFlag=DS.result;
options.pvalueMinFlag=and(pvalueMinFlag,DSI_flag);          
% fileName2=fileName;fileName2(1).name=['DSI_',fileName2(1).name];
options.fileName=fileName;
options.DSFlag=DSFlag;
colorMapforROIs_DSI(options);
%%

%%
if options.curveSave==1
    [figAll{1},ylimitSave1]=gaussianFitPlot(dataFit,dataUnfit,options);
    [figAll{2},ylimitSave2]=shadedErrorStack(dataAndError,options);
    pvalue

     saveFigs(filePath,file4,figAll,pvalue,ylimitSave1,ylimitSave2,options)    
end


  disp([filePath,'end----------------------------------------------------------'])
end
function saveFigs(filePath,file4,figAll,pvalue,ylimitSave1,ylimitSave2,input)


screenSize=get(0,'screensize');
screenWidth=screenSize(3);
result='curve';
% result='fiji_curves_12plot_overlay';
resultPath=fullfile(filePath,result);
if exist(resultPath)==7
else
    mkdir(resultPath)

end
width2=400;width1=200;

p1p2=length(figAll{1})+length(figAll{2});
stepSize=screenWidth/p1p2;
file1=fullfile(resultPath,file4);
ii=2;
p=length(figAll{ii});
     ratio=1.6357;width=width2;
     
for jj=1:p
    y1=ylimitSave2(jj,1);y2=ylimitSave2(jj,2);
     x1=ylimitSave2(jj,3);x2=ylimitSave2(jj,4);
     y1b=y1-0.05*(y2-y1);y2b=y2+0.05*(y2-y1);

    file5=strrep(file1,'.tif',['_curves',num2str(jj),'.eps']);
    figure(figAll{ii}(jj));
    set(gcf,'PaperPositionMode','auto')
    set(gcf,'position',[100+(jj-1)*stepSize   100   width   width*ratio]);
     set(gca,'TickDir','out');
     xlim([x1-5,x2]);ylim([y1b,y2b]);
      axis off;
      set(gca,'position',[0.1300    0.1100    0.7750    0.8150]);
       
       print('-r300',strrep(file5,'.eps','.png'),'-dpng')
       if input.epsSave==1
           print('-r300',file5,'-depsc','-tiff')
       end
       savefig(strrep(file5,'.eps','.fig'))
       
end

disX=100+(jj)*stepSize;
ii=1;
p=length(figAll{ii});
ratio=2.5;width=width1;
for jj=1:p
    y1=ylimitSave1(jj,1);y2=ylimitSave1(jj,2);
    y1b=y1-0.05*(y2-y1);y2b=y2+0.05*(y2-y1);

    file5=strrep(file1,'.tif',['_curvesFit',num2str(jj),'.eps']);
    figure(figAll{ii}(jj));
    set(gcf,'PaperPositionMode','auto')
    set(gcf,'position',[disX+(jj-1)*stepSize   100   width   width*ratio]);
    set(gca,'TickDir','out');
    axis off;
    set(gca,'position',[0.1300    0.1100    0.7750    0.8150]);
    xlim([-10,360]);ylim([y1b,y2b]);
      
    print('-r300',strrep(file5,'.eps','.png'),'-dpng')
       if input.epsSave==1
           print('-r300',file5,'-depsc','-tiff')
           savefig(strrep(file5,'.eps','.fig'))
       end      
%       print('-r300',strrep(file5,'.eps','.png'),'-dpng')
end
file5=strrep(file1,'.tif',['_pvalue.txt']);
save(file5,'pvalue','-ASCII');


trash=0;
end

function output=drawParameters(options)
data=options.dataSave;
if strcmp(options.prefix,'fitAngOS_noDS')||strcmp(options.prefix,'fitAng_AS_and_DS')
    data=mod(data,180);
end
options.plot_condition=1;
% options.condition=DS.result;
if options.plot_condition==1
    condition=options.condition;

else
    condition=true(size(data));
    
end
I=find(condition);
data=data(condition);    

dataSave=[I,data];
output.dataResponsive=dataSave;
if max(data)>1
    %edges=linspace(0,360,13);
    edges=0:20:360;
else
    edges=0:0.1:1;
end

if  strcmp(options.prefix,'OSI') || strcmp(options.prefix,'DSI') || strcmp(options.prefix,'pvalue') 
    data(data>=1)=1;
    %edges=linspace(0,1,10);
    edges=0:0.1:1;
elseif strcmp(options.prefix,'maxDF_F')
    try
    edges=linspace(min(data),max(data),10);
    catch me
    end
elseif ~isempty(strfind(options.prefix,'FWHM'))|| strcmp(options.prefix,'fitAngOS_noDS')|| strcmp(options.prefix,'fitAng_AS_and_DS')
    %edges=linspace(0,180,7*2);
    edges=0:10:180;
end


%FWHMsub
if strcmp(options.prefix,'OSI') || strcmp(options.prefix,'DSI') || strcmp(options.prefix,'gDSI') ...
        || strcmp(options.prefix,'gOSI')|| strcmp(options.prefix,'FWHMsub')
    conditionSplit=options.conditionSplit;
    conditionSplit(~condition)=false;
    conditionRest=~conditionSplit;
    conditionRest(~condition)=false;
    
    data=options.dataSave;
    dataSplit=data(conditionSplit);
    dataRest=data(conditionRest);
    I_Split=find(conditionSplit);
    I_Rest=find(conditionRest);
    output.dataSplit=[I_Split,dataSplit];
    output.dataRest=[I_Rest,dataRest];
%     output.dataResponsive=dataSave;
%     dataSave=[I,data];
    figure(options.fig);clf; hold on;
    fc=[0.5,0.5,0.5];
    ec=[1,0,0];
    lw=0.8*2;

     %% plot histogram of the rest of cells that are responsive;
    h1=histogram(dataRest,edges);
    %         h1=histogram(DSI,edges,'Normalization','probability');
    set(h1,'FaceColor',fc,'EdgeColor',[0.5,0.5,0.5]*0.93, 'LineWidth',lw);
%     h1.FaceColor=[0.5,0.5,0.5];
%     set(h1,'EdgeColor')
%     h1.EdgeColor=[0.5,0.5,0.5];
    %% plot histogram of cells that are OS OR DS;
    h2=histogram(dataSplit,edges);
    set(h2,'FaceColor','none','EdgeColor',ec, 'LineWidth',lw);
    x0=median(dataSplit);
    yLim=get(gca,'yLim');
%     hold on;
    plot([1,1]*x0,yLim,'r--')
    hh=text(x0,0.5*sum(yLim),num2str(x0,'%.2f'));
    hh.Color=[1 0 0];
    %%
        trash=0;
    resultPath=options.resultPath;
    fileName=options.fileName;
    file=fileName(1).name(1:11);file=[options.prefix,'_',file,'.tif'];
    set(gcf,'PaperPositionMode','auto')
    title(file,'Interpreter','none')
    print('-r300',fullfile(resultPath,strrep(file,'.tif','_hist.png')),'-dpng')
%     ylim([0,max(h2.Values)])
    clf; hold on;
        h2=histogram(dataSplit,edges);
    set(h2,'FaceColor','none','EdgeColor',ec, 'LineWidth',lw);
%     x0=median(dataSplit);
%     yLim=get(gca,'yLim');
%     hold on;
%     plot([1,1]*x0,yLim,'r--')
%     hh=text(x0,0.5*sum(yLim),num2str(x0,'%.2f'));
%     hh.Color=[1 0 0];   
    title(file,'Interpreter','none')
    print('-r300',fullfile(resultPath,strrep(file,'.tif','_hist_rescaled.png')),'-dpng')
    save(fullfile(resultPath,strrep(file,'.tif','.txt')),'dataSave','-ASCII')   
else
    
    figure(options.fig);clf;
    h1=histogram(data,edges);
    %         h1=histogram(DSI,edges,'Normalization','probability');

    h1.FaceColor=[0.5,0.5,0.5];
    h1.EdgeColor=[0.5,0.5,0.5];
    x0=median(data);
    yLim=get(gca,'yLim');
    hold on;
    plot([1,1]*x0,yLim,'b--')
    hh=text(x0,0.5*sum(yLim),num2str(x0,'%.2f'));
    hh.Color=[1 0 0];
        trash=0;
    resultPath=options.resultPath;
    fileName=options.fileName;
    file=fileName(1).name(1:11);file=[options.prefix,'_',file,'.tif'];
    set(gcf,'PaperPositionMode','auto')
    title(file,'Interpreter','none')
    print('-r300',fullfile(resultPath,strrep(file,'.tif','_hist.png')),'-dpng')
    save(fullfile(resultPath,strrep(file,'.tif','.txt')),'dataSave','-ASCII')       
end
end

function [jjAll,ylimitSave]=gaussianFitPlot(data3b,data2b,input)
% data3b, fitted datax,y; data2b,unfitted;x,y,se;
jj=100;figure(jj);jjAll=[jj];
clf;hold on;
[m,n]=size(data3b);
roiN=n-1;
x1=data3b(:,1);
x2=data2b(:,1);
s=0;
step=ones(roiN,1)*2;
step(1)=0;
t=1;
ylimitSave=[0,0];
minY2=0;
maxY2=0;
for ii=1:roiN
    if mod(ii,input.plotN)==0
        ylimitSave(t,:)=[minY2,maxY2];t=t+1;
        jj=jj+1;figure(jj);clf;hold on;
        jjAll=[jjAll;jj];
        minY2=10^5;
        maxY2=-10^5;
    end
    y1=data3b(:,ii+1);
    y2=data2b(:,2*ii);
    y2SE=data2b(:,2*ii+1);
    y2=y2-min(y2);
       %remove negtive values
        %normalize the data for fitting.
      

        maxY=max(y2);
        y2=y2/maxY;
        y2SE=y2SE/maxY;
    offSet=sum(step(1:ii));
    y2=y2-offSet;
    y1=y1-offSet;
    plot(x1,y1,'-r');
    minY2=min([minY2,min(y1)]);
    maxY2=max([maxY2,max(y1)]);
    hh2=errorbar(x2,y2,y2SE,'ob');
    yMean=max(y1);
    text(x2(end-3),yMean,['p=',num2str(input.pvalue(ii),'%04.2f')])
     text(x2(1),yMean,[num2str(ii,'%.0f')])
%     set(hh2,'MarkerSize',5)
end
 ylimitSave(t,:)=[minY2,maxY2];
end

 function [figAll,ylimitSave]=shadedErrorStack(data,input)
ylimitSave=[0,0,0,0];
minY2=0;
maxY2=0;
[m,n]=size(data);
data2=data;
roiN=(n-1)/2;% first column is x;
x=data(:,1);
s=0;
step=zeros(roiN,1);
ratio=ones(roiN,1);
for ii=2:roiN
    maxV=max(max(data(:,2:end)));
    step(ii)=maxV*1;
end
maxV2=zeros(roiN,1);
input.linepros=cell(roiN,1);
for ii=1:roiN
    dataIn=[x,data(:,2*ii:2*ii+1)*ratio(ii)];
    maxV2(ii)=max(abs(dataIn(:,2)));
    input.linepros{ii}={'-k','lineWidth',0.7};
end
% for ii=1:2:roiN
% %     input.linepros{ii}={'-r','lineWidth',0.7};
% end
% maxV2(2:2:end)=maxV2(1:2:end);
maxV2=maxV2/max(maxV2);
ratio=1./maxV2;
ratio(ratio<=2)=1;
ratio(ratio>2&ratio<=4)=2;
ratio(ratio>4&ratio<=8)=3;
ratio(ratio>8&ratio<=12)=6;
ratio(ratio>12)=10;
% ratio=ones(roiN,1);
fig=200;
figure(fig);clf;
input.fig=fig; hold on;
figAll=[fig];
 xx=1;
 scaleB=2;
%  step=step*1.2;
 t=1;
for ii=1:roiN
    if mod(ii,input.plotN)==0
        fig=fig+1;figure(fig);clf;input.fig=fig;hold on;
        figAll=[figAll;fig];
        ylimitSave(t,:)=[minY2,maxY2,x(1),x(end)];
        t=t+1;
        minY2=10^9;
        maxY2=-10^9;        
    end    
    dataIn=[x,data(:,2*ii:2*ii+1)*ratio(ii)];
    dataIn(:,2)=dataIn(:,2)-sum(step(1:ii));
    
    shadedErrorStack_sub(dataIn,input.linepros{ii});
    minY2=min([minY2,min(dataIn(:,2))]);
    maxY2=max([maxY2,max(dataIn(:,2))]);
    y0=nanmean(dataIn(:,2));
    scaleB2=scaleB*ratio(ii);
    y1=y0-scaleB2/2;
    y2=y1+scaleB2-1;
    plot([1,1]*xx,[y1,y2],input.linepros{ii}{:})
    hh=text(xx+10,0.5*y1+y2/2,[num2str(ii)]);hh.Color=[1 0 0];
end

    ylimitSave(t,:)=[minY2,maxY2,x(1),x(end)];
% ylim([-5200,200])


%  for ii=1:roiN
%      plot([1,1]*x,[50,y0+scaleB*ratio(ii)]-sum(step(1:ii)),'k-')
%  end
%  plot([1,1],[100,200],'k-')
 end

function shadedErrorStack_sub(dataIn,linepros)
x=dataIn(:,1);
p=length(x);
y=dataIn(:,2);
ySE=dataIn(:,3);
b=~isnan(ySE);
c=b(2:end)-b(1:end-1);
x1=find(c==1);x1=x1+1;
x2=find(c==-1);
x1=x1(:);
x2=x2(:);
if b(1)
    x1=[1;x1];
end
if b(end)
    x2=[x2;p];
end
xx=[x1,x2];
p1=size(xx,1);
% linepros={'-k','lineWidth',2};
for ii=1:p1
    xIn=xx(ii,1):xx(ii,2);
    shadedErrorBar(x(xIn),y(xIn),ySE(xIn),linepros);
end
trash=0;


end


function [data3,data3SE]=trialSortAvg(data,s,order)
data2=data;
[m,n]=size(data);

for ii=1:n
    data2(:,ii)=data(order(:,end),ii);
end
p2=floor(m/s);
s2=ones(m,1)*p2;

s3=(1:s).'*ones(1,p2);
s3b=s3(:);
p3=length(s3b);
s2(1:p3)=s3b(:);
input.SEflag=1;
[data3,data3SE,angleLeft2]=avgRepeatNO(data2,s2,input);
end


function b1=findDentrite(a,a1)
[m,~]=size(a);
a1=a1(:);
[m2,~]=size(a1);
b1=zeros(m2,2);
b1(:,1)=a1;
for ii=1:m2
    s=find(a(:,1)==a1(ii));
    b1(ii,2)=a(s,2);
end
end

function colorMapforROIs(input)
        map=input.map;
        avg=input.avg;  
        stdImg=input.stdImg;
        pvalueMinFlag=input.pvalueMinFlag;
        OSFlag=input.OSFlag;
        bw=input.bw;             
        I=find(OSFlag);
        fileName=input.fileName;
        resultPath=input.resultPath;
        
mapSave=map;        
[m,n]=size(avg);
map=nan(m,n);        
p2=length(I);
centroidXY=zeros(p2,2);
for ii=1:p2
    jj=I(ii);
    tmp=regionprops(bw(:,:,jj),'centroid');
    centroidXY(ii,:)=tmp.Centroid;
    map(bw(:,:,jj))= mapSave(jj);
end
%% processing avg
avg=mat2gray(avg);
avgTmp=sort(avg(:));
ratioAvg=0.99;
avgUpperLimit=avgTmp(round(ratioAvg*m*n));
avg=mat2gray(avg,[0,avgUpperLimit]);


%% 360
centroidXY(:,1)=centroidXY(:,1)-3;

[rgbImg,map2]=mapHSV(avg,map);
stdImg2=mapHSV(stdImg,map);
id=1:11;
file=[fileName(1).name(id),'_colormap.tif'];
cmbFile=['cmb_',fileName(1).name(id),'.tif'];
imwrite(uint8(rgbImg*255),fullfile(resultPath,['cmb_',file]),'tif','compression','none');
fig=25;
hf=figure(fig);clf
imshow(rgbImg,[]); 
hold on;
fontSiz=4;
for ii=1:p2
     text('position',centroidXY(ii,:),'fontsize',fontSiz,'string',num2str(I(ii)),'color',[0 1 0]) 
end
set(gcf,'PaperPositionMode','auto')
print('-r300',fullfile(resultPath,strrep(file,'.tif','Number.png')),'-dpng')
% imwrite(uint8(rgbImg*255),fullfile(resultPath,strrep(file,'.tif','Number.tif')),'tif','compression','none');
if input.epsSave==1
    print('-r300',fullfile(resultPath,strrep(file,'.tif','Number.eps')),'-depsc','-tiff')
end

%%
hf=figure(fig+1);clf
imshow(rgbImg,[]); 
hold on;
% fontSiz=8;
for ii=1:p2
     text('position',centroidXY(ii,:),'fontsize',fontSiz,'string',num2str(mapSave(I(ii)),'%.0f'),'color',[0 1 0]) 
end
set(gcf,'PaperPositionMode','auto')
print('-r300',fullfile(resultPath,strrep(file,'.tif','angle.png')),'-dpng')

if input.epsSave==1
    print('-r300',fullfile(resultPath,strrep(file,'.tif','angle.eps')),'-depsc','-tiff')
end
file=fullfile(resultPath,file);
imwrite(uint8(stdImg2*255),strrep(file,'.tif','_std.tif'),'tif','compression','none');
%%
hf=figure(fig+2);clf
imshow(stdImg2,[]); 
hold on;
for ii=1:p2
     text('position',centroidXY(ii,:),'fontsize',fontSiz,'string',num2str(I(ii)),'color',[0 1 0]) 
end
set(gcf,'PaperPositionMode','auto')
print('-r300',strrep(file,'.tif','stdNumber.png'),'-dpng')

% imwrite(uint8(stdImg2*255),strrep(file,'.tif','_stdNumber.tif'),'tif','compression','none');
%%
img=avg*0+255;
img=uint8(img);
[rgbImg,map2]=mapHSV(img,map);
[mm,nn,p]=size(rgbImg);
imgAll=zeros(mm,nn*3+4,p);
% imgAll(:,:,1)=0.5;
s=0;
s=s+1;imgAll(:,(s-1)*nn+(1:nn)+s-1,:)=rgbImg;
figure(2);clf;subplot(1,3,1);
imshow(uint8(rgbImg*256));title('360')
figure(3);subplot(1,3,1);imshow(map2)
imwrite(uint8(rgbImg*255),strrep(file,'.tif','_360Fullmap.tif'),'tif','compression','none');
%%
[rgbImg,map2]=mapHSV_180Degrees(img,map);
figure(2);
subplot(1,3,2);
imshow(uint8(rgbImg*256));title('180');
[m,n,~]=size(map2);
figure(3);subplot(1,3,2);imshow(map2(:,1:round(n/2),:))

s=s+1;imgAll(:,(s-1)*nn+(1:nn)+s,:)=rgbImg;
imwrite(uint8(rgbImg*255),strrep(file,'.tif','_180Halfmap.tif'),'tif','compression','none');
[rgbImg,map2]=mapHSV_180Degrees(avg,map);
imwrite(uint8(rgbImg*255),strrep(fullfile(resultPath,cmbFile),'.tif','_180Halfmap.tif'),'tif','compression','none');
%%
input.fullColor=1;
input.lhFlag=0;
[rgbImg,map2]=mapHSV_180Degrees(img,map,input);
imwrite(uint8(rgbImg*255),strrep(file,'.tif','_180Fullmap.tif'),'tif','compression','none');

% combined raw figure and pixel map with 180 degree full colormap
[rgbImg2,map2]=mapHSV_180Degrees(avg,map,input);
imwrite(uint8(rgbImg2*255),strrep(fullfile(resultPath,cmbFile),'.tif','_180Fullmap.tif'),'tif','compression','none');


s=s+1;imgAll(:,(s-1)*nn+(1:nn)+s,:)=rgbImg;
figure(2);
subplot(1,3,3);
imshow(uint8(rgbImg*256));title('180');
[m,n]=size(map2);
figure(3);subplot(1,3,3);imshow(map2);
file5=strrep(file,'.tif','_blank.eps');
figure(2);
print('-r300',strrep(file5,'.eps','.png'),'-dpng')
set(gcf,'PaperPositionMode','auto')
print('-r300',file5,'-depsc','-tiff')
figure(3);
file5=strrep(file,'.tif','_clcbar.eps');
print('-r300',strrep(file5,'.eps','.png'),'-dpng');
set(gcf,'PaperPositionMode','auto')
print('-r300',file5,'-depsc','-tiff')
imwrite(uint8(imgAll*255),strrep(file,'.tif','_blk_img.tif'),'tif','compression','none');
end

 function colorMapforROIs_DSI(input)
        map=input.map;
        avg=input.avg;  
        stdImg=input.stdImg;
        pvalueMinFlag=input.pvalueMinFlag;
        bw=input.bw;     
        DSFlag=input.DSFlag;
        %I=find(pvalueMinFlag);
        I=find(DSFlag);
        fileName=input.fileName;
        resultPath=input.resultPath;
        
mapSave=map;        
[m,n]=size(avg);
map=nan(m,n);        
p2=length(I);
centroidXY=zeros(p2,2);
for ii=1:p2
    jj=I(ii);
    tmp=regionprops(bw(:,:,jj),'centroid');
    centroidXY(ii,:)=tmp.Centroid;
    map(bw(:,:,jj))= mapSave(jj);
end


%%
centroidXY(:,1)=centroidXY(:,1)-3;

[rgbImg,map2]=mapHSV(avg,map);
stdImg2=mapHSV(stdImg,map);
id=1:11;
file=['DSI_',fileName(1).name(id),'_colormap.tif'];
imwrite(uint8(rgbImg*255),fullfile(resultPath,file),'tif','compression','none');
fig=25;
hf=figure(fig);clf
imshow(rgbImg,[]); 
hold on;
fontSiz=4;
for ii=1:p2
     text('position',centroidXY(ii,:),'fontsize',fontSiz,'string',num2str(I(ii)),'color',[0 1 0]) 
end
set(gcf,'PaperPositionMode','auto')
print('-r300',fullfile(resultPath,strrep(file,'.tif','Number.png')),'-dpng')
TRASH=0;
% imwrite(uint8(rgbImg*255),fullfile(resultPath,strrep(file,'.tif','Number.tif')),'tif','compression','none');
% print('-r300',fullfile(resultPath,strrep(file,'.tif','Number.eps')),'-depsc','-tiff')
%%
hf=figure(fig+1);clf
imshow(rgbImg,[]); 
hold on;
% fontSiz=8;
for ii=1:p2
     text('position',centroidXY(ii,:),'fontsize',fontSiz,'string',num2str(mapSave(I(ii)),'%.0f'),'color',[0 1 0]) 
end
set(gcf,'PaperPositionMode','auto')
print('-r300',fullfile(resultPath,strrep(file,'.tif','angle.png')),'-dpng')
% print('-r300',fullfile(resultPath,strrep(file,'.tif','angle.eps')),'-depsc','-tiff')

file=fullfile(resultPath,file);
imwrite(uint8(stdImg2*255),strrep(file,'.tif','_std.tif'),'tif','compression','none');
%%
hf=figure(fig+2);clf
imshow(stdImg2,[]); 
hold on;
for ii=1:p2
     text('position',centroidXY(ii,:),'fontsize',fontSiz,'string',num2str(I(ii)),'color',[0 1 0]) 
end
set(gcf,'PaperPositionMode','auto')
print('-r300',strrep(file,'.tif','stdNumber.png'),'-dpng')

% imwrite(uint8(stdImg2*255),strrep(file,'.tif','_stdNumber.tif'),'tif','compression','none');
%
img=avg*0+255;
img=uint8(img);
[rgbImg,map2]=mapHSV(img,map);
[mm,nn,p]=size(rgbImg);
imgAll=zeros(mm,nn*3+4,p);
% imgAll(:,:,1)=0.5;
s=0;
s=s+1;imgAll(:,(s-1)*nn+(1:nn)+s-1,:)=rgbImg;
figure(2);clf;subplot(1,3,1);
imshow(uint8(rgbImg*256));title('360')
figure(3);subplot(1,3,1);imshow(map2)
imwrite(uint8(rgbImg*255),strrep(file,'.tif','_360Fullmap.tif'),'tif','compression','none');
 end
%%
% [rgbImg,map2]=mapHSV_180Degrees(img,map);
% figure(2);
% subplot(1,3,2);
% imshow(uint8(rgbImg*256));title('180');
% [m,n,~]=size(map2);
% figure(3);subplot(1,3,2);imshow(map2(:,1:round(n/2),:))
% 
% s=s+1;imgAll(:,(s-1)*nn+(1:nn)+s,:)=rgbImg;
% imwrite(uint8(rgbImg*255),strrep(file,'.tif','_180Halfmap.tif'),'tif','compression','none');
% %%
% input.fullColor=1;
% input.lhFlag=0;
% [rgbImg,map2]=mapHSV_180Degrees(img,map,input);
% imwrite(uint8(rgbImg*255),strrep(file,'.tif','_180Fullmap.tif'),'tif','compression','none');
% s=s+1;imgAll(:,(s-1)*nn+(1:nn)+s,:)=rgbImg;
% figure(2);
% subplot(1,3,3);
% imshow(uint8(rgbImg*256));title('180');
% [m,n]=size(map2);
% figure(3);subplot(1,3,3);imshow(map2);
% file5=strrep(file,'.tif','_blank.eps');
% figure(2);
% print('-r300',strrep(file5,'.eps','.png'),'-dpng')
% set(gcf,'PaperPositionMode','auto')
% print('-r300',file5,'-depsc','-tiff')
% figure(3);
% file5=strrep(file,'.tif','_clcbar.eps');
% print('-r300',strrep(file5,'.eps','.png'),'-dpng');
% set(gcf,'PaperPositionMode','auto')
% print('-r300',file5,'-depsc','-tiff')
% imwrite(uint8(imgAll*255),strrep(file,'.tif','_blk_img.tif'),'tif','compression','none');