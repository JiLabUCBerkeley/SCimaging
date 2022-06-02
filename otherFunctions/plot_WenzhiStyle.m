function plot_WenzhiStyle(data,order,options)
    
    %[order,qt]=StimulationSequence(filePath,fileName,p);
if nargin==3
    options.qt=[13,12,10];options.greenDot=1:3;
    options.trims=[3,0];
end
trialNO=qt(3);
angleNO=qt(2);
framePerSti2=length(order)/angleNO/trialNO;
order2=reshape(order/(360/angleNO)+1,angleNO,trialNO);
order2=round(order2);
[DFF,f0s,baseline]=calculatedDff0(data,1);
sorted_intensity=reshape(data,qt(1),qt(2),qt(3));
sorted_dff0=reshape(DFF,qt(1),qt(2),qt(3));
for irep=1:1:qt(3)
    sorted_intensity(:,order2(:,irep),irep,:)=sorted_intensity(:,:,irep,:);
    sorted_dff0(:,order2(:,irep),irep,:)=sorted_dff0(:,:,irep,:);
end
%%
figure(1), clf;
set(gcf,'PaperType','usletter','Unit','inches');
set(gcf,'PaperPosition',[0.5 0.5 7.5 10]);
ax(1)=subplot(20,4,1:12);
plot_intensity(ax(1),data,f0s);   
hold on;
qt=options.qt;
qt23=qt(2)*qt(3);
qt1=qt(1);
yTrash=zeros(qt1,qt23);
yTrash(options.greenDot,:)=1;
yTrash(end-options.trims(1):end,:)=2;
 yTrash(options.trims(1)+1,:)=3;
yTrash2=yTrash(:);
xTrash=1:length(yTrash2);
xTrash=xTrash(:);
I_green=find(yTrash2==1);
I_blue=find(yTrash2==2);
I_red=find(yTrash2==3);
set(ax(1),'NextPlot','add'); 
plot(ax(1),xTrash(I_green),Intensity(I_green,iROI),'.g')
plot(ax(1),xTrash(I_red),Intensity(I_red,iROI),'.r')
xlabel(['green:',num2str(options.greenDot(1)),'-',num2str(options.greenDot(end)),';Red:',num2str(options.trims(1)+1)])
%% dff
input=options;
        ax(2)=subplot(20,4,17:28);
        plot_dff(ax(2),DFF,f0s);
                ax(3)=subplot(20,4,[37 38 41 42 45 46 49 50 53 54]);
        plot_tuningcurve(ax(3),sorted_dff0(:,:,:,1),qt,[5 0],'normal',10);
        
        ax(4)=subplot(20,4,[39 40 43 44 47 48 51 52 55 56]);hold on;
        %%
    deleteN=options.trims;
    order=deleteMN(order,deleteN,framePerSti2);
    [Y,I]=sort(order(:,1));
    order(:,2:4)=[(1:length(I)).',Y(:)*1,I(:)]; 
    framePerSti=length(order)/angleNO/trialNO;
    qt=[framePerSti,angleNO,trialNO];
    options.angleNO=angleNO;
    options.trialNO=trialNO;
    options.framePerSti=framePerSti;
    options.qt=qt;
    y3=deleteMN(data,deleteN,framePerSti2);
    pvalue=pvalueGet(y3(order(:,end)),options); 
    output=GaussianFit(y3,order);
    ii=1;p1=1;
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
    dataFit(:,1+ii)=output.z2(:);
    dataUnfit(:,2*ii:2*ii+1)=[output.y1(:),output.y1SE(:)];  
    dataAndError(:,2*ii:2*ii+1)=[data3b(:),data3SEb(:)];

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
        title(ax(4), ['Pvalue=' num2str(pvalue) ', \theta=' num2str( output.theta) ', SSE =' num2str(output.SSE) ', RSQ =' num2str(utput.RSQ)]);
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
         polar([xdata_interp+90 xdata_interp(1)+90]/180*pi, [y1 y1(1)], 'b');
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
         
         is_goodROI='No';
         if (maxydata(iROI)>=goodROI_threshold), is_goodROI='Yes'; end
         text(0,0.95,sprintf('This ROI is a good ROI(Max respons>=%d%%): %s\n max response=%.2f',goodROI_threshold,is_goodROI,maxydata(iROI)),'Parent',ax(6));
         isgoodfitting='No';
         if (OS_fitting_sse(iROI)<=SSE_Threshold)&&...
            (OS_fitting_rsq(iROI)>=RSquared_Threshold)&&...
            (maxydata(iROI)>=goodROI_threshold)
            isgoodfitting='Yes';
         end
         text(0,0.8,sprintf('A good fitting(fitting SSE<%.2f AND RSQ>%.2f): %s',SSE_Threshold,RSquared_Threshold,isgoodfitting),'Parent',ax(6));
         
         isOrientationTuned='No';
         if (Pvalue_anova(iROI)<Pvalue_anova_threshold)&&strcmp(isgoodfitting,'Yes'),isOrientationTuned='Yes'; end
         text(0,0.6,sprintf('This ROI is OS(anova test pValue<%.3f): %s',Pvalue_anova_threshold,isOrientationTuned),'Parent',ax(6));
         osstr=sprintf('OS: Theta=%.2f, FWHM=%.2f,\n gOSI=%.2f, OSI=%.2f, DSI=%.2f',estimates(2),estimates(3)*FWHM_const,gOSI(iROI),OSI(iROI),DSI(iROI));
         if strcmp(isOrientationTuned,'No')
           osstr=sprintf('OS: Theta=%.2f, FWHM=%.2f,\n gOSI=%.2f, OSI=%.2f, DSI=%.2f',NaN,NaN,NaN,NaN,NaN);
         else
             OT(iROI)=1;
         end
         text(0,0.45,osstr,'Parent',ax(6));
         
         isDS='No';
         if (DSI(iROI)>=0.5)&&strcmp(isOrientationTuned,'Yes'), isDS='Yes'; end
         text(0,0.3,sprintf('This ROI is DS(DSI>=0.5): %s',isDS),'Parent',ax(6));
         DSstr=sprintf('DS: Theta=%.2f, gDirection=%.2f, FWHM=%.2f,\n gDSI=%.2f, DSI=%.2f',estimates(2),gDirection(iROI)*180/pi,estimates(3)*FWHM_const,gDSI(iROI),DSI(iROI));
         if strcmp(isDS,'No')
             DSstr=sprintf('DS: Theta=%.2f, gDirection=%.2f, FWHM=%.2f,\n gDSI=%.2f, DSI=%.2f',NaN,NaN,NaN,NaN,NaN);
         else
            DS(iROI)=1;
         end
         text(0,0.15,DSstr,'Parent',ax(6));
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
