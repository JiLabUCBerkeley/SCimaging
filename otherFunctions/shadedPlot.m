function shadedPlot(dff0,pAll,dim,options)
    % dff0 is m by n dimension. n is the ROI quantity; pAll indicates how
    % to organize the data, e.g., pAll=[p1,p2, p3, p4]; Typically, p1 means
    % framePerstimulus; p2 means the quantity of temproal frequency; p3
    % means the quantity of spatial frequency, p4 means the trialNO; dim
    % indicates the dimension to do averaging. 
    [m,n]=size(dff0);
[avg,SE]=avgStd(dff0,pAll,dim);
dimLeft=1:length(pAll);
dimLeft(dim)=[];
pLeft=pAll(dimLeft);
data2b=[(1:size(avg,1)).',avg,SE];
data2b(:,2:2:end)=avg;
data2b(:,3:2:end)=SE;
v.pLeft=pLeft;
v.dim=dim;
v.pAll=pAll;
v.dimLeft=dimLeft;
options.v=v;
options=plotCurves_SE_v2(data2b,options);
figAll=options.figAll;
x=data2b(:,1);
%% add information on it;
% if pLeft(1)==1
% x2=reshape(x,pLeft(1),prod(pLeft(2:end)));
% 
% % x2Mean=mean(x2);
% x2Mean=0.5*x2(1,:)+0.5*mean(x2);
% x2End=x2(end,:);
% orderAll=options.orderAll;
% orderAll2=avgStd(orderAll,pAll,dim);
% orderAll3=orderAll2(:,dimLeft);
% dimInfo=zeros(length(pLeft),length(x2Mean));
% if length(pLeft)>2
%     x3=reshape(x,prod(pLeft(1:2)),prod(pLeft(3:end)));
%     x3End=x3(end,:);
% end
% for ii=1:length(pLeft)
%     x=orderAll3(:,ii);
%     x2=reshape(x,pLeft(1),prod(pLeft(2:end)));
%     x2=x2(1,:)/min(x2(1,:));
%     c=unique(x2);
%     c=sort(c);
%     for jj=1:length(c)
%         x2(x2==c(jj))=jj;
%     end
%     dimInfo(ii,:)=x2;
% end
% stringLeft=options.dimension(dimLeft);
% for ii=1:length(figAll)
%     h=get(figAll(ii),'CurrentAxes');
%     yMax=h.YLim(2);
%     
%     for jj=1:length(x2Mean)
%         printText=cell(length(stringLeft)-1,1);
%         for ss=2:length(stringLeft)
%             printText{ss-1}=[stringLeft{ss},num2str(dimInfo(ss,jj),'%.0f')];
%             hText=text(h,x2Mean(jj),yMax-(ss-2)*200,printText{ss-1},'color',[1,0,1]);
%             hText.FontSize=10;
% %             text(h,x2Mean(jj),yMax,printText{ss-1},'color',[1,0,1])
%         end
%         
%         
%         plot(h,[1,1]*x2End(jj),h.YLim,'--','color',[1,1,1]*0.7,'LineWidth',0.5)
%         
%     end
%     
%     if length(pLeft)>2
%        for jj=1:length(x3End)
%            plot(h,[1,1]*x3End(jj),h.YLim,'-','color',[1,0,0]*1,'LineWidth',0.5)
%        end
% %         x3End=x3(end,:);
%     end
% end
% end