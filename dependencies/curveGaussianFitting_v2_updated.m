function output=curveGaussianFitting_v2(x,y,SE,output)
% minOptions=2; % 1 is wenzhi's way: move the whole curve to positive; 2 force negative values to zero

p=length(y);
ydata=y;output.minY=min(ydata);
ydata = double(ydata);

%normalize the data for fitting.
ii=1;

maxY=max(ydata);
SE=SE/maxY;
ydata=ydata/maxY;
ydata=[ydata(end); ydata; ydata(1)];
SE=[SE(end); SE; SE(1)];
%to fit normal with a double guassian function
%Initiate fitting parameters
x=x(:);
xdata=[x(1)-x(2);x;x(2)-x(1)+x(end)];
xdata=xdata.';
xdata_interp=0:1:360;
dc = min(ydata);
[amp1, amp1_idx] = max(ydata);
%sigma1 = 15;
%sigma2 = 15;
sigma1 = 50;
sigma2 = 50;
theta = xdata(amp1_idx);
amp2_idx = amp1_idx-6;
if amp2_idx < 1
   amp2_idx = amp2_idx + 8;
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
resultPath = 'E:\2p_DATA\Suite2p_results\2018_05_02\ANM409927\combined_s13_s14_s15\session15';
%save(fullfile(resultPath, 'estimates'))
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

if estimates(2)<0, estimates(2) = 360+estimates(2); end
estimates(2) = rem(estimates(2),360);
theta=estimates(2);
output.theta=theta;

output.z2=z_interp;

output.x2=xdata_interp(:);
%% DSI flag
z = ori_gauss(estimates, [estimates(2) rem(estimates(2)+180,360) rem(estimates(2)+90,360) rem(estimates(2)+270,360)]);
if output.negativeMethod==1
    z(z<0)=0;
elseif output.negativeMethod==2
    if min(z)<0
        z=z-min(z);
    end
end
OSI_ii=(z(1)+z(2)-z(3)-z(4))/sum(z); 
DSI_ii=(z(1)-z(2))/(z(1)+z(2));
output.DSI_ii=DSI_ii;
output.OSI_ii=OSI_ii;
output.FWHM=estimates(3)*2.3548;
output.SSE=SSE;
output.RSQ=RSQ;

function [estimates, SSE, resnorm, residual,exitflag, lambda, jacobian] = fitcurve_OriGauss(xdata, ydata, start_point, LB, UB)
options = optimset('Display','off');
[estimates,resnorm, residual,exitflag, output, lambda, jacobian] = lsqcurvefit(@ori_gauss, start_point, xdata, ydata, LB, UB,options);
z = ori_gauss(estimates, xdata);
SSE = sum((z-ydata).*(z-ydata));
