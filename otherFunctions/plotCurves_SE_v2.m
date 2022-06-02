function options=plotCurves_SE_v2(data2b,options)
% data3b, fitted datax,y; data2b,unfitted;x,y,se;
if nargin==1
    options.pupilFlag=0;
    options.plotN=12;
    options.scaleY=50;
    options.scaleX=20;
    options.stepScale=1;
    options.filePath=cd;
    options.file='curve';
    options.curveSave=1;
    options.epsSave=0;
end
data=data2b(:,2:2:end);
t=data2b(:,1);
dataSE=data2b(:,3:2:end);
[m,n]=size(data);

plotN=options.plotN;
figN=floor(n/plotN);
idNTmp=cell(figN,1);
for ii=1:figN
    idNTmp{ii}=(1:plotN)+(ii-1)*plotN;
end
if mod(n,plotN)>0
    figN=figN+1;
    idN=cell(figN,1);
    idN(1:end-1)=idNTmp;
    if isempty(idNTmp)
        idN{figN}=1:n;
    else
        idN{figN}=idNTmp{end}(end)+1:n;
    end
    
else
    idN=idNTmp;
end
%%
data2=cell(figN,1);
data2SE=cell(figN,1);

if options.pupilFlag==1
    if isfield(options,'pupil')
        pupilFlag=1;
        pupil=options.pupil;
    else
        pupilFlag=0;
    end    
else
    pupilFlag=0;
end
%% sort correlation
if pupilFlag==1
    R2=corrcoef([pupil(:),data]);
    R=R2(2:end,1);
    [~,I]=sort(R,'descend');
%     [R_sort,I]=sort(abs(R),'decendent');
%     [~,I]=sort(abs(R),'descend');
    R_sort=R(I);
    data=data(:,I);
    R_save=corrcoef([pupil(:),data]);
    I_save=[I,R_sort];
    idN2=idN;
    R3=idN;
    b=0;
    for ii=1:length(idN2)
        p1=length(idN2{ii});
        b=b+p1;
        idN2{ii}=I(b-p1+1:b).';
        R3{ii}=R_sort(b-p1+1:b).';
    end
else
    idN2=idN;
    

end

for ii=1:figN
    data2{ii}=data(:,idN{ii});
    data2SE{ii}=dataSE(:,idN{ii});
    if pupilFlag==1
        data2SE{ii}=[pupil,data2SE{ii}];
    end
end

figAll=1:figN;figAll=figAll+99;
x=t;xx=0.5;
     scaleY=options.scaleY;
     scaleX=options.scaleX; 
     ylimitSave=zeros(figN,4);
for ii=1:figN
    figure(figAll(ii));clf; hold on;
    data3=data2{ii};
    data3SE=data2SE{ii};
    [m,roiN]=size(data3);
    
    step=zeros(roiN,1);
    ratio=ones(roiN,1);
    maxV=max(max(abs(data3(:,2:end))));

    maxV2=max(data3,[],1);
    maxV2=maxV2/max(maxV2);
    ratio=1./maxV2;
    ratio(ratio<=2)=1;
    ratio(ratio>2&ratio<=4)=2;
    ratio(ratio>4&ratio<=8)=3;
    ratio(ratio>8&ratio<=12)=6;
    ratio(ratio>12)=10;  
    ratio=ratio*0+1;
    options.linepros=cell(roiN,1);
    options.linepros{1}={'color',[1,1,1]*0.4,'lineWidth',0.7};
    options.linepros(2:end)=options.linepros(1);
   
    if pupilFlag==1
        
        options.linepros{1}={'-r','lineWidth',1.2};
        ratio(1)=3;
    end
     for jj=2:roiN
%          minY_pre=min(data3(:,jj)*ratio(jj));
%          maxY_cur=max(data3(:,jj)*ratio(jj));
%         step(jj)=(maxY_cur-minY_pre)*options.stepScale;
        step(jj)=maxV*options.stepScale;
    end   
    for jj=1:roiN
        data3(:,jj)=data3(:,jj)*ratio(jj)-sum(step(1:jj));
        data3SE(:,jj)=data3SE(:,jj)*ratio(jj)*1;
        shadedErrorBar(x,data3(:,jj),data3SE(:,jj),'-k');
%          data3SE(:,jj)=data3SE(:,jj)*ratio(jj)-sum(step(1:jj));
%         plot(x,data3(:,jj),options.linepros{jj}{:})
%         errorbar(x,data3(:,jj),data3SE(:,jj),'-b');
        y0=nanmean(data3(:,jj));
        scaleY2=scaleY*ratio(jj);
        y1=y0-scaleY2/2;
        y2=y1+scaleY2-1;
        plot([1,1]*xx,[y1,y2],'r','lineWidth',0.7)
        if jj==1
%             plot(-[0,scaleX]+x(end),(0.5*y1+y2/2)*[1,1],'-g','lineWidth',0.7)
%             hh=text(x(end)-scaleX,0.5*y1+y2/2+0.5*(y2-y1),[num2str(scaleX),' s']);hh.Color=[1 0 0];
            minY2_1=min(data3(:,jj));
            maxY2_1=max(data3(:,jj));
        elseif jj==2
            hh=text(xx+10,0.5*(minY2_1+max(data3(:,jj))),[num2str(scaleY),'%\DeltaF/F']);hh.Color=[0 0 0];

        end
         jj2=jj;
    if pupilFlag==1
        jj2=jj-1;
        if jj==1
            hh=text(xx+10,0.5*y1+y2/2,'pupil');hh.Color=[1 0 1];
        else
            plot(x(options.f0s{idN2{ii}(jj2)}),data3(options.f0s{idN2{ii}(jj2)},jj),'b.')
            hh=text(xx+10,0.5*y1+y2/2,[num2str(idN2{ii}(jj2)),',R=',num2str(R3{ii}(jj2),'%.2f')]);hh.Color=[1 0 1];
        end
    else
%         plot(x(options.f0s{idN2{ii}(jj2)}),data3(options.f0s{idN2{ii}(jj2)},jj),'b.')
        hh=text(xx+10,0.5*y1+y2/2,[num2str(idN2{ii}(jj2))]);hh.Color=[1 0 1];
        

    end   
    hh.FontSize=14;
        
        
    end
    minY2=min(data3(:));
    maxY2=max(data3(:));
   ylimitSave(ii,:)=[minY2,maxY2,min([x(1),xx]),x(end)];
     
%     maxV2=zeros(roiN,1);
    
end 

screenSize=get(0,'screensize');
screenWidth=screenSize(3);
result='curve';
% result='fiji_curves_12plot_overlay';
resultPath=fullfile(options.filePath,result);
if exist(resultPath)==7
else
    mkdir(resultPath)

end
options.resultPath=resultPath;
width2=400;width1=150;
options.figAll=figAll;
p1p2=length(figAll);
stepSize=screenWidth/p1p2;
file=options.file;
file=['avg_',file];
file1=fullfile(resultPath,file);
if pupilFlag==1
    R=im2uint8(mat2gray(double(R_save),[-1,1]));
    if size(R,1)<20
        R=imresize(R,10,'nearest');
    elseif size(R,1)<40
        R=imresize(R,5,'nearest');
    elseif size(R,1)<100
        R=imresize(R,4,'nearest');
    elseif size(R,1)<300
         R=imresize(R,2,'nearest');
    end
    
    figure(3);clf;imshow(R);colormap(jet(256));colorbar;
    colorSelected=jet(256);
    imgColor=ind2rgb(R,colorSelected);
    imwrite(imgColor,strrep(file1,'.tif','cor.tif'),'tif','compression','none')
    R=ones(20,1)*linspace(-1,1,128);R=im2uint8(mat2gray(double(R),[-1,1]));
    imgColor=ind2rgb(R,colorSelected);
    imwrite(imgColor,strrep(file1,'.tif','map.tif'),'tif','compression','none')    
    save(strrep(file1,'.tif','roi_cor.txt'),'I_save','-ASCII')
end
ii=1;
p=length(figAll);
     ratio=1.6357;width=width2;
%%
% figAll=options.figAll;
x=data2b(:,1);
v=options.v;

pLeft=v.pLeft;
dim=v.dim;
pAll=v.pAll;
dimLeft=v.dimLeft;
%% add information on it;
% if pLeft(1)==1
x2=reshape(x,pLeft(1),prod(pLeft(2:end)));

% x2Mean=mean(x2);
x2Mean=0.5*x2(1,:)+0.5*mean(x2);
x2End=x2(end,:);
orderAll=options.orderAll;
orderAll2=avgStd(orderAll,pAll,dim);
orderAll3=orderAll2(:,dimLeft);
dimInfo=zeros(length(pLeft),length(x2Mean));
if length(pLeft)>2
    x3=reshape(x,prod(pLeft(1:2)),prod(pLeft(3:end)));
    x3End=x3(end,:);
end
for ii=1:length(pLeft)
    x=orderAll3(:,ii);
    x2=reshape(x,pLeft(1),prod(pLeft(2:end)));
    if min(x2(1,:))==0
        x2=x2(1,:);
    else
        x2=x2(1,:)/min(x2(1,:));
    end
    
    c=unique(x2);
    c=sort(c);
    for jj=1:length(c)
        x2(x2==c(jj))=jj;
    end
    dimInfo(ii,:)=x2;
end
stringLeft=options.dimension(dimLeft);
for ii=1:length(figAll)
    h=get(figAll(ii),'CurrentAxes');figure(figAll(ii));
    yMax=h.YLim(2);
    
    for jj=1:length(x2Mean)
        printText=cell(length(stringLeft)-1,1);
        for ss=2:length(stringLeft)
            printText{ss-1}=[stringLeft{ss},num2str(dimInfo(ss,jj),'%.0f')];
            hText=text(x2Mean(jj),yMax-(ss-2)*200,printText{ss-1},'color',[1,0,1]);
            hText.FontSize=10;
%             text(h,x2Mean(jj),yMax,printText{ss-1},'color',[1,0,1])
        end
        
        
        plot(h,[1,1]*x2End(jj),h.YLim,'--','color',[1,1,1]*0.7,'LineWidth',0.5)
        
    end
    
    if length(pLeft)>2
       for jj=1:length(x3End)
           plot(h,[1,1]*x3End(jj),h.YLim,'-','color',[1,0,0]*1,'LineWidth',0.5)
       end
%         x3End=x3(end,:);
    end
end
%%

for jj=1:p
    y1=ylimitSave(jj,1);y2=ylimitSave(jj,2);
     x1=ylimitSave(jj,3);x2=ylimitSave(jj,4);
     y1b=y1-0.05*(y2-y1);y2b=y2+0.05*(y2-y1);

    file5=strrep(file1,'.tif',['_curves',num2str(jj),'.eps']);
    figure(figAll(jj));
    set(gcf,'PaperPositionMode','auto')
    set(gcf,'position',[5+(jj-1)*stepSize   100   width*2   width*ratio]);
     set(gca,'TickDir','out');
     xlim([x1,x2]);ylim([y1b,y2b]);
%      xLine=options.xLine;
%      for ss=1:length(xLine)
%          plot([1,1]*xLine(ss),[y1b,y2b],'--r')
%      end
      axis off;
      set(gca,'position',[0.1300    0.1100    0.7750    0.8150]);
       if options.curveSave==1
           print('-r300',strrep(file5,'.eps','.png'),'-dpng')
       end
       
       if options.epsSave==1
           print('-r300',file5,'-depsc','-tiff')
       end
%        savefig(strrep(file5,'.eps','.fig'))
       
end