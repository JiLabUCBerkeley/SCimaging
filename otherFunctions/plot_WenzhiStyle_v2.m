function plot_WenzhiStyle_v2(data,dff0,f0s,order,qt,options)
    
    %[order,qt]=StimulationSequence(filePath,fileName,p);
% if nargin==3
%     options.qt=[13,12,10];options.greenDot=1:3;
%     options.deleteN=[3,0];
% end
options.greenDot=1:options.deleteN(1);
if isempty(options.greenDot)
    options.greenDot=1;
end
angleNO=qt(3);
trialNO=qt(2);

framePerSti2=length(order)/angleNO/trialNO;
qt=[framePerSti2,angleNO,trialNO];qt_raw=qt;
order2=reshape(order(1:framePerSti2:end,1)/(360/angleNO)+1,angleNO,trialNO);
order2=round(order2);
% [dff0,f0s,baseline]=calculatedDff0(data,1);
sorted_intensity=reshape(data,qt(1),qt(2),qt(3));
sorted_dff0=reshape(dff0,qt(1),qt(2),qt(3));
for irep=1:1:qt(3)
    sorted_intensity(:,order2(:,irep),irep,:)=sorted_intensity(:,:,irep,:);
    sorted_dff0(:,order2(:,irep),irep,:)=sorted_dff0(:,:,irep,:);
end
%%
p1=1;
warning('off','stats:statrobustfit:IterationLimit');
%mod
deleteN=options.deleteN;
nanROI=options.nanROI;
% angleNO=qt(3);
% trialNO=qt(2);
framePerSti2=length(order)/angleNO/trialNO;
order=deleteMN(order,deleteN,framePerSti2);
[Y,I]=sort(order(:,1));
order(:,2:4)=[(1:length(I)).',Y(:)*1,I(:)]; 
framePerSti=length(order)/angleNO/trialNO;
qt=[framePerSti,angleNO,trialNO];
options.angleNO=angleNO;
options.trialNO=trialNO;
options.framePerSti=framePerSti;
options.qt=qt;
% pvalueMin=0.01;
map=nan(p1,1);
pvalue=ones(p1,1);
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
             pvalue(ii)=pvalueGet(y3(order(:,end)),options); 
             %%
             options.negativeMethod=0;
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


             %%
              options.SEflag=1;
              [data3,data3SE,angleLeft2]=avgTrials(y3,order(:,1),framePerSti,angleNO,trialNO,options);
              data3b=insertNan(data3,framePerSti,nanRatio);
              data3SEb=insertNan(data3SE,framePerSti,nanRatio);
              %%

             if ii==1
                 xdata_interp=output.x1;
                 x2=output.x2;
                 dataFit=zeros(length(x2),p1+1);
                 dataFit(:,1)=x2(:);
                 dataUnfit=zeros(length(xdata_interp),1+2*p1);
                 dataFit(:,1+ii)=output.z2(:);
                 dataUnfit(:,2*ii:2*ii+1)=[output.y1(:),output.y1SE(:)];
                 dataUnfit(:,1)=xdata_interp(:);
                 
                 dataAndError=zeros(length(data3b),1+2*p1);
                 dataAndError(:,1)=(1:length(data3b)).';
                 dataAndError(:,2*ii:2*ii+1)=[data3b(:),data3SEb(:)];
             end
                 dataFit(:,1+ii)=output.z2(:);
                 dataUnfit(:,2*ii:2*ii+1)=[output.y1(:),output.y1SE(:)];  
                 dataAndError(:,2*ii:2*ii+1)=[data3b(:),data3SEb(:)];
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
        angle_rawTmp=find(xTrash==min(xTrash));
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

% %%
% result='curve';
% % result='fiji_curves_12plot_overlay';
% resultPath1=fullfile(filePath,result);
% if exist(resultPath1)==7
% else
%     mkdir(resultPath1)
% 
% end
% result='map';
% resultPath=fullfile(resultPath1,result);
% if exist(resultPath)
% else
%     mkdir(resultPath)
% end
% 
% result='params';
% resultPath2=fullfile(resultPath1,result);
% if exist(resultPath2)
% else
%     mkdir(resultPath2)
% end


% options.resultPath=resultPath2;
% options.fileName=fileName;
options.pvalue=pvalue;
options.DSI=DSI;
options.OSI=OSI;
options.gDSI=gDSI;
%%
figure(1), clf;
set(gcf,'PaperType','usletter','Unit','inches');
set(gcf,'PaperPosition',[0.5 0.5 7.5 10]);
ax(1)=subplot(20,4,1:12);
plot_intensity(ax(1),data,f0s);   
hold on;
qt=options.qt;
qt23=qt(2)*qt(3);
qt1=framePerSti2;
yTrash=zeros(qt1,qt23);
yTrash(options.greenDot,:)=1;
yTrash(end-options.deleteN(1):end,:)=2;
 yTrash(options.deleteN(1)+1,:)=3;
yTrash2=yTrash(:);
xTrash=1:length(yTrash2);
xTrash=xTrash(:);
I_green=find(yTrash2==1);
I_blue=find(yTrash2==2);
I_red=find(yTrash2==3);
set(ax(1),'NextPlot','add'); 

iROI=1;
plot(ax(1),xTrash(I_green),data(I_green,iROI),'.g')
plot(ax(1),xTrash(I_red),data(I_red,iROI),'.r')
xlabel(['green:',num2str(options.greenDot(1)),'-',num2str(options.greenDot(end)),';Red:',num2str(options.deleteN(1)+1)])
%% dff
input=options;
        ax(2)=subplot(20,4,17:28);
        plot_dff(ax(2),dff0,f0s);
                ax(3)=subplot(20,4,[37 38 41 42 45 46 49 50 53 54]);
        plot_tuningcurve(ax(3),sorted_dff0(:,:,:,1),qt_raw,[5 0],'normal',10);
        
        ax(4)=subplot(20,4,[39 40 43 44 47 48 51 52 55 56]);hold on;
        %%
%     deleteN=options.deleteN;
%     order=deleteMN(order,deleteN,framePerSti2);
%     [Y,I]=sort(order(:,1));
%     order(:,2:4)=[(1:length(I)).',Y(:)*1,I(:)]; 
%     framePerSti=length(order)/angleNO/trialNO;
%     qt=[framePerSti,angleNO,trialNO];
%     options.angleNO=angleNO;
%     options.trialNO=trialNO;
%     options.framePerSti=framePerSti;
%     options.qt=qt;
%     y3=deleteMN(data,deleteN,framePerSti2);
%     pvalue=pvalueGet(y3(order(:,end)),options); 
%     output=GaussianFit(y3,order);
    ii=1;p1=1;
%     xdata_interp=output.x1;
%     x2=output.x2;
%     dataFit=zeros(length(x2),p1+1);
%     dataFit(:,1)=x2(:);
%     dataUnfit=zeros(length(xdata_interp),1+2*p1);
%     dataFit(:,1+ii)=output.z2(:);
%     dataUnfit(:,2*ii:2*ii+1)=[output.y1(:),output.y1SE(:)];
%     dataUnfit(:,1)=xdata_interp(:);
%     dataAndError=zeros(length(data3b),1+2*p1);
%     dataAndError(:,1)=(1:length(data3b)).';
%     dataAndError(:,2*ii:2*ii+1)=[data3b(:),data3SEb(:)];
%     dataFit(:,1+ii)=output.z2(:);
%     dataUnfit(:,2*ii:2*ii+1)=[output.y1(:),output.y1SE(:)];  
%     dataAndError(:,2*ii:2*ii+1)=[data3b(:),data3SEb(:)];

    xdata_interp=dataFit(:,1);
    x2=dataUnfit(:,1);
    s=0;
    step=ones(1,1)*2;
    step(1)=0;
    t=1;
    ylimitSave=[0,0];
    minY2=0;
    maxY2=0;
    ii=1;
    y1=dataFit(:,ii+1);
    y2=dataUnfit(:,2*ii);
    y2SE=dataUnfit(:,2*ii+1);
    y2=y2-min(y2);
    %remove negtive values
    %normalize the data for fitting.


    maxY=max(y2);
    y2=y2/maxY;
    y2SE=y2SE/maxY;
    offSet=0;
    y2=y2-offSet;
    y1=y1-offSet;
    plot(xdata_interp,y1,'-r');
    minY2=min([minY2,min(y1)]);
    maxY2=max([maxY2,max(y1)]);
    hh2=errorbar(x2,y2,y2SE,'ob','Parent',ax(4));
    maxydata=max(y2);iROI=1;

 
%         hold off;
%         xlim(ax(4), [xdata_interp(1) xdata_interp(end)]);
        set(ax(4), 'XTick', 0:360/(input.qt(2)):360);
        set(ax(4), 'TickDir', 'out');
        xlabel(ax(4), 'Degree(\circ)');
        ylabel(ax(4), 'Normalized \DeltaF/F');
        xlim(ax(4), [-5 365]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        title(ax(4), ['Pvalue=' num2str(pvalue) ', \theta=' num2str( output.theta) ', SSE =' num2str(output.SSE) ', RSQ =' num2str(output.RSQ)]);
        yl=get(ax(4),'YLim');
        xdata=(0:(input.qt(2)-1))*360/(input.qt(2));
         ax(5)=subplot(20,4,[63 64 67 68 71 72 75 76 79 80]);
         dori = diff(xdata(1:2));
         xdata = rem(xdata+90,360)/180*pi;
         mw = max(max(y2/max(y2)), max(y1));
         r = circ_r(xdata,(y2/maxydata(1))',dori) * mw;
         phi = circ_mean(xdata,(y2/maxydata(1))');
         %     xdata_interp = xdata_interp/180*pi;
         hold on
         zm = r*exp(i*phi');
         
         % draw a unit circle
         zz = exp(i*linspace(0, 2*pi, 101)) * mw;
         plot(real(zz),imag(zz),'-', 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5])
         plot(real(zz/2),imag(zz/2),'-','LineWidth', 0.5, 'Color', [0.5 0.5 0.5])
         plot(real(0/2),imag(0/2),'-','LineWidth', 0.5, 'Color', [0.5 0.5 0.5])
         plot([-mw mw], [0 0], '-', [0 0], [-mw mw], '-', 'Color', [0.5 0.5 0.5])
         %
         polar([xdata xdata(1)], [(y2/maxydata(iROI))' y2(1)/maxydata(iROI)], 'ok')
         polar([xdata_interp+90; xdata_interp(1)+90]/180*pi, [y1; y1(1)], 'b');
         hold off
         
         formatSubplot(gca,'ax','square','box','off','lim',[-mw-0.1 mw+0.1 -mw-0.1 mw+0.1]);
         set(gca,'xtick',[]);
         set(gca,'ytick',[]);
         text(-max(real(zz))/20,max(real(zz))*1.1,'Up');
         text(-max(real(zz))/9,-max(real(zz)*1.1),'Down');
         text(max(real(zz))*1.1,max(real(zz))/5,'Posterior','Rotation',270);
         text(-max(real(zz))*1.1,-max(real(zz))/6,'Anterior','Rotation',90);
         axis off;
         
         ax(6)=subplot(20,4,[61 62 65 66 69 70 73 74 77 78]);
         set(ax(6),'XLim',[0 1],'Ylim',[0 1]);
         axis off;
          text(0,0.6,['fittedAngle=',num2str(map,'%.0f'), '; pvalue=',num2str(pvalue,'%.4f'),';'],'Parent',ax(6));
          text(0,0.3,[OS.output.flagHead{1},num2str(OS.output.OS_flagAll(1)),';  ',DS.output.flagHead{1},num2str(DS.output.DS_flagAll(1)) ],'Parent',ax(6));
%           text(0,0.8,sprintf('A good fitting(fitting SSE<%.2f AND RSQ>%.2f): %s',SSE_Threshold,RSquared_Threshold,isgoodfitting),'Parent',ax(6));
%           text(0,0.6,sprintf('This ROI is OS(anova test pValue<%.3f): %s',Pvalue_anova_threshold,isOrientationTuned),'Parent',ax(6));
%          is_goodROI='No';
%          if (maxydata(iROI)>=goodROI_threshold), is_goodROI='Yes'; end
%          text(0,0.95,sprintf('This ROI is a good ROI(Max respons>=%d%%): %s\n max response=%.2f',goodROI_threshold,is_goodROI,maxydata(iROI)),'Parent',ax(6));
%          text(0,0.45,osstr,'Parent',ax(6));
%          text(0,0.3,osstr,'Parent',ax(6));
%          text(0,0.15,osstr,'Parent',ax(6));
%          isgoodfitting='No';
%          if (OS_fitting_sse(iROI)<=SSE_Threshold)&&...
%             (OS_fitting_rsq(iROI)>=RSquared_Threshold)&&...
%             (maxydata(iROI)>=goodROI_threshold)
%             isgoodfitting='Yes';
%          end
%          text(0,0.8,sprintf('A good fitting(fitting SSE<%.2f AND RSQ>%.2f): %s',SSE_Threshold,RSquared_Threshold,isgoodfitting),'Parent',ax(6));
%          
%          isOrientationTuned='No';
%          if (Pvalue_anova(iROI)<Pvalue_anova_threshold)&&strcmp(isgoodfitting,'Yes'),isOrientationTuned='Yes'; end
%          text(0,0.6,sprintf('This ROI is OS(anova test pValue<%.3f): %s',Pvalue_anova_threshold,isOrientationTuned),'Parent',ax(6));
%          osstr=sprintf('OS: Theta=%.2f, FWHM=%.2f,\n gOSI=%.2f, OSI=%.2f, DSI=%.2f',estimates(2),estimates(3)*FWHM_const,gOSI(iROI),OSI(iROI),DSI(iROI));
%          if strcmp(isOrientationTuned,'No')
%            osstr=sprintf('OS: Theta=%.2f, FWHM=%.2f,\n gOSI=%.2f, OSI=%.2f, DSI=%.2f',NaN,NaN,NaN,NaN,NaN);
%          else
%              OT(iROI)=1;
%          end
%          text(0,0.45,osstr,'Parent',ax(6));
%          
%          isDS='No';
%          if (DSI(iROI)>=0.5)&&strcmp(isOrientationTuned,'Yes'), isDS='Yes'; end
%          text(0,0.3,sprintf('This ROI is DS(DSI>=0.5): %s',isDS),'Parent',ax(6));
%          DSstr=sprintf('DS: Theta=%.2f, gDirection=%.2f, FWHM=%.2f,\n gDSI=%.2f, DSI=%.2f',estimates(2),gDirection(iROI)*180/pi,estimates(3)*FWHM_const,gDSI(iROI),DSI(iROI));
%          if strcmp(isDS,'No')
%              DSstr=sprintf('DS: Theta=%.2f, gDirection=%.2f, FWHM=%.2f,\n gDSI=%.2f, DSI=%.2f',NaN,NaN,NaN,NaN,NaN);
%          else
%             DS(iROI)=1;
%          end
%          text(0,0.15,DSstr,'Parent',ax(6));
function formatSubplot(handle,varargin)

args.fs = 7;
args.xl = [];
args.yl = [];
args.box = [];
args.ax = [];
args.lim = [];
args.tt = [];
args.xt = [];
args.yt =[];
args = parseVarArgs(args,varargin{:});

set(handle,'fontsize',args.fs)
if ~isempty(args.xl)
  xlabel(handle,args.xl)
end
if ~isempty(args.yl)
  ylabel(handle,args.yl)
end
if ~isempty(args.tt)
  title(handle,args.tt)
end
if ~isempty(args.box)
  set(handle,'box',args.box)
end
if ~isempty(args.ax)
  axis(handle,args.ax)
end
if ~isempty(args.lim)
  axis(handle,args.lim)
end
if ~isempty(args.yt)
  set(handle,'ytick',args.yt)
end
if ~isempty(args.xt)
  set(handle,'xtick',args.xt)
end
function [output_args]=plot_intensity(ahandle,intensity,f0s)
%Inputs:    ahandle
%           intensity and indics of data points selected as F0 as baseline
%           of dff calculation.
%created by Wenzhi Sun, 08.17.2015
if nargin<3, f0s=[]; end
if ishandle(ahandle)
    cla(ahandle);
    Imax = max(intensity);
    Imin = min(intensity);
    xdata = 1:length(intensity);
    hold on
    plot(ahandle,xdata, intensity, 'k', 'LineWidth', .5);
%     plot(ahandle,f0s, intensity(f0s), 'm.', 'MarkerSize', 6, 'MarkerFaceColor', 'none');
    hold off;
    box('off');
    set(ahandle, 'TickDir', 'out');
    set(ahandle, 'XLim', [xdata(1)-25 xdata(end)+25]);
    set(ahandle, 'YLim', [Imin*0.95 Imax+Imin/20]);
    set(ahandle, 'FontName', 'arial', 'FontSize', 10);
    
    ylabel(ahandle,'Fluo. Intensity (a.u.)','FontName','arial','FontAngle','normal','FontSize',12);
    xlabel(ahandle,'Frame #','FontName','arial','FontSize',12);
end        

function [output_args]=plot_dff(ahandle,dff,f0s)
%Inputs:    ahandle
%           dff and indics of f0s 
%created by Wenzhi Sun, 08.17.2015
if nargin<3, f0s=[]; end
if ishandle(ahandle)
    cla(ahandle);
    dmax = max(dff);
    dmin = min(dff);
    xdata = 1:length(dff);
    hold on
    plot(ahandle,xdata, dff, 'k', 'LineWidth', .5);
    plot(ahandle,f0s, dff(f0s), 'bx', 'MarkerSize', 6, 'MarkerFaceColor', 'none');
    hold off;
    box('off');
    set(ahandle, 'TickDir', 'out');
    set(ahandle, 'XLim', [xdata(1)-25 xdata(end)+25]);
    set(ahandle, 'YLim', [dmin*0.95 dmax+dmin/20]);
    set(ahandle, 'FontName', 'arial', 'FontSize', 10);
    
    ylabel(ahandle,'\DeltaF/F0','FontName','arial','FontAngle','italic','FontSize',12);
    xlabel(ahandle,'Frame #','FontName','arial','FontSize',12);  
end
function output = plot_tuningcurve(ahandle,sorted_data,ndimens,trims,plot_title,nblanks)

if nargin<6, nblanks=10;end
if nargin<5, plot_title='normal'; end
if nargin<4, trims=[0 0]; end

output = 0;

data_plot = nan(ndimens(1)+nblanks, ndimens(2), ndimens(3));
data_plot(1:ndimens(1),:,:) = sorted_data;
data_mean_rep = reshape(data_plot, ndimens(2)*(ndimens(1)+nblanks), ndimens(3));
data_se_rep = nanstd(data_mean_rep,0,2)/sqrt(ndimens(3));
data_mean_rep = nanmean(data_mean_rep,2);

data_mean_dir =squeeze(nanmean(sorted_data(trims(1)+1:end-trims(2),:,:),1));
ydata_mean=nanmean(data_mean_dir,2);
yerror_mean=nanstd(data_mean_dir,0,2)/sqrt(ndimens(3));
if min(ydata_mean)<0
    data_mean_dir=data_mean_dir-min(yerror_mean);
    yerror_mean=yerror_mean-min(yerror_mean);
    yerror_mean=yerror_mean-min(yerror_mean);
end

hold on
xdata = 1:1:ndimens(2)*(ndimens(1)+nblanks);
errorbar(xdata(1:end-nblanks+1), data_mean_rep(1:end-nblanks+1), data_se_rep(1:end-nblanks+1),'Color',[0.7 0.7 0.7],'LineWidth',0.5, 'Parent', ahandle);

xdata_mean=(0:ndimens(2)-1)*(ndimens(1)+nblanks)+round(ndimens(1)/2);
errorbar(ahandle,xdata_mean,ydata_mean,yerror_mean,'LineStyle','none','Marker', 'o','MarkerSize',10,'LineWidth',2, 'Parent', ahandle);
hold off;

box('off');
set(ahandle, 'TickDir', 'out');
set(ahandle, 'XLim', [xdata(1)-25 xdata(end)+25],'XTick', xdata_mean, 'XTickLabel',(1:ndimens(2))-1);
set(ahandle, 'FontName', 'arial', 'FontSize', 10);

ylabel(ahandle,'\DeltaF/F0','FontName','arial','FontAngle','italic','FontSize',12);
xlabel(ahandle,'Test #','FontName','arial','FontSize',12);
function r = circ_r(alpha, w, d, dim)
% r = circ_r(alpha, w, d)
%   Computes mean resultant vector length for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [w		number of incidences in case of binned angle data]
%     [d    spacing of bin centers for binned data, if supplied 
%           correction factor is used to correct for bias in 
%           estimation of r, in radians (!)]
%     [dim  compute along this dimension, default is 1]
%
%     If dim argument is specified, all other optional arguments can be
%     left empty: circ_r(alpha, [], [], dim)
%
%   Output:
%     r		mean resultant length
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N.I. Fisher
%   Topics in circular statistics, S.R. Jammalamadaka et al. 
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if nargin < 4
  dim = 1;
end

if nargin < 2 || isempty(w) 
  % if no specific weighting has been specified
  % assume no binning has taken place
	w = ones(size(alpha));
else
  if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1) 
    error('Input dimensions do not match');
  end 
end

if nargin < 3 || isempty(d)
  % per default do not apply correct for binned data
  d = 0;
end

% compute weighted sum of cos and sin of angles
r = sum(w.*exp(1i*alpha),dim);

% obtain length 
r = abs(r)./sum(w,dim);

% for data with known spacing, apply correction factor to correct for bias
% in the estimation of r (see Zar, p. 601, equ. 26.16)
if d ~= 0
  c = d/2/sin(d/2);
  r = c*r;
end

function [mu ul ll] = circ_mean(alpha, w, dim)
%
% mu = circ_mean(alpha, w)
%   Computes the mean direction for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [w		weightings in case of binned angle data]
%     [dim  compute along this dimension, default is 1]
%
%     If dim argument is specified, all other optional arguments can be
%     left empty: circ_mean(alpha, [], dim)
%
%   Output:
%     mu		mean direction
%     ul    upper 95% confidence limit
%     ll    lower 95% confidence limit 
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N. I. Fisher
%   Topics in circular statistics, S. R. Jammalamadaka et al. 
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if nargin < 3
  dim = 1;
end

if nargin < 2 || isempty(w)
  % if no specific weighting has been specified
  % assume no binning has taken place
	w = ones(size(alpha));
else
  if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1) 
    error('Input dimensions do not match');
  end 
end

% compute weighted sum of cos and sin of angles
r = sum(w.*exp(1i*alpha),dim);

% obtain mean by
mu = angle(r);

% confidence limits if desired
if nargout > 1
  t = circ_confmean(alpha,0.05,w,[],dim);
  ul = mu + t;
  ll = mu - t;
end

function t = circ_confmean(alpha, xi, w, d, dim)
%
% t = circ_mean(alpha, xi, w, d, dim)
%   Computes the confidence limits on the mean for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [xi   (1-xi)-confidence limits are computed, default 0.05]
%     [w		number of incidences in case of binned angle data]
%     [d    spacing of bin centers for binned data, if supplied 
%           correction factor is used to correct for bias in 
%           estimation of r, in radians (!)]
%     [dim  compute along this dimension, default is 1]
%
%   Output:
%     t     mean +- d yields upper/lower (1-xi)% confidence limit
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N. I. Fisher
%   Topics in circular statistics, S. R. Jammalamadaka et al. 
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if nargin < 5
  dim = 1;
end

if nargin < 4 || isempty(d)
  % per default do not apply correct for binned data
  d = 0;
end

if nargin < 3 || isempty(w)
  % if no specific weighting has been specified
  % assume no binning has taken place
	w = ones(size(alpha));
else
  if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1) 
    error('Input dimensions do not match');
  end 
end

% set confidence limit size to default
if nargin < 2 || isempty(xi)
  xi = 0.05;
end

% compute ingredients for conf. lim.
r = circ_r(alpha,w,d,dim);
n = sum(w,dim);
R = n.*r;
c2 = chi2inv((1-xi),1);

% check for resultant vector length and select appropriate formula
t = zeros(size(r));

for i = 1:numel(r)
  if r(i) < .9 && r(i) > sqrt(c2/2/n(i))
    t(i) = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));  % equ. 26.24
  elseif r(i) >= .9
    t(i) = sqrt(n(i)^2-(n(i)^2-R(i)^2)*exp(c2/n(i)));      % equ. 26.25
  else 
    t(i) = NaN;
    warning('Requirements for confidence levels not met.');
  end
end

% apply final transform
t = acos(t./R);
function params = parseVarArgs(params,varargin)
% Parse variable input arguments supplied in name/value format.
%
%    params = parseVarArgs(params,'property1',value1,'property2',value2) sets
%    the fields propertyX in p to valueX.
%
%    params = parseVarArgs(params,varargin{:},'strict') only sets the field
%    names already present in params. All others are ignored.
%
% AE 2007-06-01

if isempty(varargin)
    return
end

% check if correct number of inputs
if mod(length(varargin),2)
    if ~strcmp(varargin{end},'strict')
        err.message = 'Name and value input arguments must come in pairs.';
        err.identifier = 'parseVarArgs:wrongInputFormat';
        error(err)
    else
        % in 'strict' case, remove all fields that are not already in params
        fields = fieldnames(params);
        ndx = find(~ismember(varargin(1:2:end-1),fields));
        varargin([2*ndx-1 2*ndx end]) = [];
    end
end

% parse arguments
for i = 1:2:length(varargin)
    if ischar(varargin{i})
        params.(varargin{i}) = varargin{i+1};
    else
        err.message = 'Name and value input arguments must come in pairs.';
        err.identifier = 'parseVarArgs:wrongInputFormat';
        error(err)
    end
end