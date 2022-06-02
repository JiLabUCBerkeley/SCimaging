function t1=temporalFilter(t1,filterMethod)
%% filter shifts
% 
% filterMethod.name=1; % 0: none; 1: gaussian; 2: median; 3: Wiener
% filterMethod.sizNo=7;
% filterMethod.sigmaNo=5; % no used for median & wiener filters
% curve=IOS_time_gui_filter(curve,filterMethod,'vertical');


name=filterMethod.name;

if name==0
else 
    sizeNumber=filterMethod.sizeNumber;
    if name==1;
        sigmaNumber=filterMethod.sigmaNumber;    
        w=fspecial('gaussian',[sizeNumber,1],sigmaNumber);
        t1=imfilter(t1,w,'replicate');
    elseif name==2
        t1=medianFilterAdapt(t1,sizeNumber);
    elseif name==3
                    t1=wiener2(t1,[sizeNumber,1]);
    elseif name==4
        sigmaNumber=filterMethod.sigmaNumber;   
        t1=edge(t1,'canny',[.05,.4],sigmaNumber);
    end
    
end
% 
% 
% [m1,n1]=size(t1);
% if nargin==2
%     if n1>m1
%         t1=t1.';
%     end
%     if name==0
%         %imgAvg=imgAvg;
%     elseif name==1
%         w=fspecial('gaussian',[sizeNumber,1],sigmaNumber);
%         t1=imfilter(t1,w,'replicate');
%     elseif name==2
%         t1=IOS_time_medianFilterAdapt(t1,sizeNumber);
%     elseif name==3
%             t1=wiener2(t1,[sizeNumber,1]);
%     elseif name==4
%         t1=edge(t1,'canny',[.05,.4],sigmaNumber);
% 
%     end
% 
%     if n1>m1
%         t1=t1.';
%     end
% elseif nargin==3
% 
%     if name==0
%         %imgAvg=imgAvg;
%     elseif name==1
%         w=fspecial('gaussian',[sizeNumber,1],sigmaNumber);
%         t1=imfilter(t1,w,'replicate');
%     elseif name==2
%         t1=medianFilterAdapt(t1,sizeNumber);
%     elseif name==3
%             t1=wiener2(t1,[sizeNumber,1]);
%     elseif name==4
%         t1=edge(t1,'canny',[.05,.4],sigmaNumber);
% 
%     end
% 
%     
% end

function x=medianFilterAdapt(x,n)
[m,n2]=size(x);
y=zeros(m+2*n,n2);
for ii=1:n2
    y(n+1:n+m,ii)=x(:,ii);
    y(1:n,ii)=x(1,ii);
    y(n+m+1:end,ii)=x(end,ii);
%     y(:,ii)=medfilt1(y(:,ii),n);
    y(:,ii)=medfilt2(y(:,ii),[n,1]);
    x(:,ii)=y(n+1:n+m,ii);   
end



