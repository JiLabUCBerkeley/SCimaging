function z = ori_gauss(params, xdata)
%use this double gaussian function to fit orientation tuning curve
%     amp1= params(1);
%     theta = params(2);
%     sigma1 = params(3);
%     amp2 = params(4);
%     sigma2 = params(5);
%     dc = params(6);
% created by Wenzhi, 04/10/2013

%different sigma
% z =params(1)*exp(-((-360*(sign(xdata-params(2)-180)>0)+360*(sign(xdata-params(2)+180)<0)+xdata-params(2)).^2)/params(3)/params(3)/2)...
%     +params(4)*exp(-((360*(sign(xdata-params(2))<0)+xdata-params(2)-180).^2)/params(5)/params(5)/2)...
%     +params(6);

%same sigma
z =params(1)*exp(-((-360*(sign(xdata-params(2)-180)>0)+360*(sign(xdata-params(2)+180)<0)+xdata-params(2)).^2)/params(3)/params(3)/2)...
    +params(4)*exp(-((360*(sign(xdata-params(2))<0)+xdata-params(2)-180).^2)/params(3)/params(3)/2)...
    +params(6);

return

 