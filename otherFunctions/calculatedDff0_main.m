function [dff0,f0,baseline]=calculatedDff0_main(thisIntensity,method,input)
% method 1, mod; most frequent mod
%method 2 fixed: fixPersti has two inputs; fixPersti(1) is the static frame number per stimulus;
% fixPersti(2) is the moving frame number per stimulus;

% method 3: fix is a vector, defing the frame numbers as baseline; example:
% fix=[1,10,11,12], then average the frame # 1, 10, 11, 12 as baseline

% method 4:  use the minimum 20% values as baseline;
[m,n]=size(thisIntensity);
dff0=zeros(m,n);
f0=cell(n,1);
baseline=zeros(n,1);
if nargin==2
    for ii=1:n
        [dff0(:,ii),f0{ii},baseline(ii)]=calculatedDff0(thisIntensity(:,ii),method);
    end    
elseif nargin==3
    for ii=1:n
        [dff0(:,ii),f0{ii},baseline(ii)]=calculatedDff0(thisIntensity(:,ii),method,input);
    end    
    
end
