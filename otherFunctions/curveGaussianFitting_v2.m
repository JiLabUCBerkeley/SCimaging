function output=curveGaussianFitting_v2(x,y,SE,output)
% minOptions=2; % 1 is wenzhi's way: move the whole curve to positive; 2 fource negative values as zero
%x=[0:30:330]
p=length(y);
% output.minOptions=minOptions;
        ydata=y;output.minY=min(ydata);
        ydata = double(ydata);
%         if min(ydata)<0
%             if minOptions==1
%                 minFlag=1;
%                 minY=min(ydata);
%                 ydata=ydata-min(ydata);
%                 
%             elseif minOptions==2
%                 minFlag=1;
%                 ydataMinFlag=ydata(1:end-2)<0;
%                 if sum(ydataMinFlag)>=(p-2)/2
%                     output.minOptions=1;
%                     minY=min(ydata);
%                     ydata=ydata-min(ydata);
%                 else
%                     ydata(ydata<0)=0;  
%                     output.minOptions=2;
%                 end
%                             
%             end
% 
%         else
%             minFlag=0;
%             output.minOptions=0;
%        end %remove negtive values
        %normalize the data for fitting.
        ii=1;

        maxY=max(ydata);
        SE=SE/maxY;
        ydata=ydata/maxY;
        ydata=[ydata(end); ydata; ydata(1)];
        SE=[SE(end); SE; SE(1)];
        %to fit normal with a double guassian function
        %Initiate fitting parameters
%         xdata=((0:numel(ydata)-1)-1)*(360/input.qt(2));
        x=x(:);
        xdata=[x(1)-x(2);x;x(2)-x(1)+x(end)];
        xdata=xdata.';
        xdata_interp=0:1:360;
        dc = min(ydata);
        [amp1, amp1_idx] = max(ydata);
        sigma1 = 70;%was 15, changed to 70 for SC on sept 8, 2017
        sigma2 = 70;%was 15, changed to 70 for SC on sept 8, 2017
        theta = xdata(amp1_idx);
        %%20180601 comment line 53 to line 56, add line 57 to 61
%         amp2_idx = amp1_idx-4;%amp2_idx = amp1_idx-6;
%         if amp2_idx < 1
%            amp2_idx = amp2_idx + 8;% amp2_idx = amp2_idx + 12;
%         end
        if amp1_idx<=length(ydata)/2
            amp2_idx=amp1_idx+length(ydata)/2;
        else
            amp2_idx=amp1_idx-length(ydata)/2;
        end
        amp2 = ydata(amp2_idx);
        
        start_point = [amp1 theta sigma1 amp2 sigma2 dc];
        LB=output.LB;
        UB=output.UB;
%         LB = [1 -30 5 0.00001 5 0];
%         UB = [1.5 360 180 1.5 180 .2];
        UB(6) = nanmean(ydata);
        
        %fitting
        [estimates, SSE, resnorm, residual, exitflag, lambda, jacobian] = fitcurve_OriGauss(xdata, ydata', start_point, LB, UB);
        
        %if first peak < second peak, exchange the theta's and redo
        %fitting ,which can make sue theta is always coressonding to
        %the biger peak at preferred orientation.
        if estimates(1)<estimates(4)
            tmpamp = estimates(1);
            estimates(1) = estimates(4);
            estimates(4) = tmpamp;
            estimates(2) = rem(estimates(2)+180, 360);
            [estimates, SSE, resnorm, residual, exitflag, lambda, jacobian] = fitcurve_OriGauss(xdata, ydata', estimates, LB, UB);
        end
        SSresid = sum(residual.^2);
        SStotal = (length(ydata)-1) * var(ydata);
        RSQ = 1 - SSresid/SStotal;
        
        z = ori_gauss(estimates, xdata);
        z_interp = ori_gauss(estimates, xdata_interp);
        
%         %Goodness of Fitting
%         OS_fitting_exitflag(ii) = exitflag;
%         %SSE
%         OS_fitting_sse(ii) = SSE;
%         %coefficient of determination, r_squred
%         SSresid = sum(residual.^2);
%         SStotal = (length(ydata)-1) * var(ydata);
%         OS_fitting_rsq(ii) = 1 - SSresid/SStotal;
%         
%         %fittings
%         data_fitted_interp(:,ii) = z_interp;
        if estimates(2)<0, estimates(2) = 360+estimates(2); end
        estimates(2) = rem(estimates(2),360);
        theta=estimates(2);
        output.theta=theta;
        
        output.z2=z_interp;
%         if minFlag==1 && normFlag==0
%             output.z2=z_interp(:)*maxY;
%             output.z2=output.z2+minY;
%         else
%             output.z2=z_interp;
%         end
     output.x2=xdata_interp(:);
%% DSI flag
z = ori_gauss(estimates, [estimates(2) rem(estimates(2)+180,360) rem(estimates(2)+90,360) rem(estimates(2)+270,360)]);
if output.negativeMethod==1
    z(z<0)=0;
elseif output.negativeMethod==2
    if min(z)<0
        z=z-min(z);
    end
else
    
end
OSI_ii=(z(1)+z(2)-z(3)-z(4))/sum(z); 
DSI_ii=(z(1)-z(2))/(z(1)+z(2));
output.DSI_ii=DSI_ii;
output.OSI_ii=OSI_ii;
output.FWHM=estimates(3)*2.3548;
output.SSE=SSE;
output.RSQ=RSQ;
% output.minOptions=minOptions;

% output.minY=min(ydata);
%         OS_fitting_AMP1(ii) = estimates(1);
%         OS_fitting_Theta(ii) = estimates(2);
%         OS_fitting_Sigma(ii) = estimates(3);
%         OS_fitting_AMP2(ii) = estimates(4);
%         OS_fitting_DC(ii) = estimates(6);
        
%         z = ori_gauss(estimates, [estimates(2) rem(estimates(2)+180,360) rem(estimates(2)+90,360) rem(estimates(2)+270,360)]);
%         OSI(ii) = (z(1)+z(2)-z(3)-z(4))/sum(z);
%         DSI(ii)=(z(1)-z(2))/(z(1)+z(2));
%         %             % % plot for visualization
%         %             if SSE<=SSE_Threshold && OS_fitting_rsq(ii)>=RSquared_Threshold
%         figure(2), clf
%         ytheta=[0 OS_fitting_AMP1(ii)+OS_fitting_DC(ii)];
%         xdirection=[gDirection(ii) gDirection(ii)]*180/pi;
%         xtheta=[OS_fitting_Theta(ii) OS_fitting_Theta(ii)];
%         hold on;
%         errorbar(xdata(2:end-1), ydata(2:end-1), SE(2:end-1), 'ob');
%         plot(xdata_interp, z_interp, '-c');
%         plot(xtheta,ytheta,'r');
%         plot(xdirection,ytheta,'g');
%         hold off;
%         xlim(gca, [xdata_interp(1) xdata_interp(end)]);
%         set(gca, 'XTick', 0:(360/input.qt(2)):360);
%         set(gca, 'TickDir', 'out');
%         xlabel(gca, 'Degree(\circ)');
%         ylabel(gca, 'Normalized \DeltaF/F');
%         xlim(gca, [-5 365]);
%         title(gca, ['Orentation normal, Pvalue=' num2str(Pvalue_anova(ii)) ', \theta=' num2str( estimates(2)) ', fiting SSE =' num2str(SSE) ', RSQ =' num2str(OS_fitting_rsq(ii)) ', OSI=' num2str(OSI(ii))]);
%         yl=get(gca,'YLim');
%         projname_print=strrep(projname_print,'_','-');
%         text(5,yl(1)+0.05,[projname_print ':' num2str(ii)]);
function [estimates, SSE, resnorm, residual,exitflag, lambda, jacobian] = fitcurve_OriGauss(xdata, ydata, start_point, LB, UB)

% opt=optimset('MaxFunEvals', 90000000, 'MaxIter', 900000000, 'TolFun',  1.000000e-08, 'Display','on');
% [estimates,resnorm, residual,exitflag, output, lambda, jacobian] = lsqcurvefit(@ori_gauss, start_point, xdata, ydata, LB, UB, opt);
% z = ori_gauss(estimates, xdata);
% SSE = sum((z-ydata).*(z-ydata));
% return;

% changed by Rongwen Lu on 2016.02.16 to speed up
% opt=optimset('MaxFunEvals', 90000000, 'MaxIter', 900000000, 'TolFun',  1.000000e-08, 'Display','on');
options = optimset('Display','off');
[estimates,resnorm, residual,exitflag, output, lambda, jacobian] = lsqcurvefit(@ori_gauss, start_point, xdata, ydata, LB, UB,options);
z = ori_gauss(estimates, xdata);
SSE = sum((z-ydata).*(z-ydata));
