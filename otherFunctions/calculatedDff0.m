function [dff0,f0,baseline,baseline1]=calculatedDff0(thisIntensity,method,input)
% method 1, mod; most frequent mod
%method 2 fixed: fixPersti has two inputs; fixPersti(1) is the static frame number per stimulus;
% fixPersti(2) is the moving frame number per stimulus;

% method 3: fix is a vector, defing the frame numbers as baseline; example:
% fix=[1,10,11,12], then average the frame # 1, 10, 11, 12 as baseline

% method 4:  use the minimum 20% values as baseline;
if method==1
    [N, X] = hist(thisIntensity, 50);
    f0_mode = find(N==max(N));
    if length(f0_mode) > 1, f0_mode = f0_mode(1); end
    half_width = (X(end)-X(1))/(length(X)-1)/2;
    f0=find(thisIntensity>X(f0_mode)-half_width&thisIntensity<X(f0_mode)+half_width);
    baseline = nanmean(thisIntensity(f0));
    dff0 = (thisIntensity - baseline)/baseline*100;
    baseline1=baseline;
elseif method==2
    % fixPersti has two inputs; fixPersti(1) is the static frame number per stimulus;
    % fixPersti(2) is the moving frame number per stimulus;
        fixPersti=input.fixPersti;
        blankII=fixPersti(1);
        stiII=fixPersti(2);
        p=length(thisIntensity);
        x1=zeros(p,1);
        for ss=1:blankII
            x1(ss:stiII:end)=1;
            
        end
        f0=find(x1==1);
        baseline = nanmean(thisIntensity(f0));
        dff0 = (thisIntensity - baseline)/baseline*100;   
            baseline1=baseline;
elseif method==3
    % fix is a vector, defing the frame numbers as baseline; example:
    % fix=[1,10,11,12], then average the frame # 1, 10, 11, 12 as baseline
        fix=input.fix;

        baseline = nanmean(thisIntensity(fix));
        dff0 = (thisIntensity - baseline)/baseline*100;   
            baseline1=baseline;
elseif method==4
    % use the minimum 20% values as baseline;
    [Y,I] = sort(thisIntensity,'ascend');
    I=I(:);
    ratio2=.2;%%change the percentile 10% or 20%
    x0=I(1:round(ratio2*length(thisIntensity)));
    f0=x0;
    baseline = nanmean(thisIntensity(f0));
    dff0 = (thisIntensity - baseline)/baseline*100;  
        baseline1=baseline;
elseif method==5
    baseline=input.baseFixed;
%     thisIntensity=(thisIntensity-baseline);
    dff0 = (thisIntensity - baseline)/baseline*100;   
    [Y,I] = sort(abs(dff0),'ascend');
    I=I(:);
    ratio2=.1;%%change the percentile 10% or 20%
    x0=I(1:round(ratio2*length(thisIntensity)));
    f0=x0;
        baseline1=baseline;
elseif method==6
    %%moving average as baseline;
    baseline1 = movmean(thisIntensity,30);
    dff0 = (thisIntensity - baseline1)./baseline1*100;
    [N, X] = hist(thisIntensity, 50);%get f0
    f0_mode = find(N==max(N));
    if length(f0_mode) > 1, f0_mode = f0_mode(1); end
    half_width = (X(end)-X(1))/(length(X)-1)/2;
    f0=find(thisIntensity>X(f0_mode)-half_width&thisIntensity<X(f0_mode)+half_width);
    baseline = nanmean(thisIntensity(f0));
elseif method==7
    %%moving average as baseline;
    [m,n]=size(thisIntensity);
    fs=m/40;
    dff0_onefs=zeros(fs,40);
    mm=reshape(thisIntensity,fs,40);
    for ii=1:40
        aa=mm(:,ii);
        y=movmean(aa,movemeanN);
      %calculate dff0, remove background (baseline)
        cc = (aa - y)./y*100;      
        dff0_onefs(:,ii)=cc;
    end 
    dff0 = reshape(dff0_onefs,m,1);
    [N, X] = hist(thisIntensity, 50);%get f0
    f0_mode = find(N==max(N));
    if length(f0_mode) > 1, f0_mode = f0_mode(1); end
    half_width = (X(end)-X(1))/(length(X)-1)/2;
    f0=find(thisIntensity>X(f0_mode)-half_width&thisIntensity<X(f0_mode)+half_width);
    baseline = nanmean(thisIntensity(f0));  
    
    
end