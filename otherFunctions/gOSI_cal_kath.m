function [gOSI,exception]=gOSI_cal_kath(ydata,datapoints,options)
if nargin==2
    options.negativeMethod=2;
end
if ~isfield(options,'negativeMethod')
    options.negativeMethod=0;
end
%Inputs: ydata a vector of responseses and how many datapoints in this
%measurement. numel(ydata)==datapoints
%created by Wenzhi Sun, July 27 2015
    if isvector(ydata)
        if eq(numel(ydata),datapoints)
            if iscolumn(ydata), ydata=ydata'; end
            theta=0:(2*pi/datapoints):(2*pi-0.1);
            if options.negativeMethod==2
                if min(ydata)<0, ydata=ydata-min(ydata); end %remove the negtive value
            elseif options.negativeMethod==1
                ydata(ydata<0)=0;
            end
            
            gOSI=abs(sum(ydata.*exp(1i*2*theta))/sum(ydata));
            exception=[];
        else
            gOSI=nan;
            exception=MException(['ws_IMAGEBOX::' 'THE DATA POINTS OF INPUT DATA DOSN''T MATCH THE INPUT DATAPOINTS!'],'Iputs Not Match');
        end
    else
        gOSI=nan;
        exception=MException(['ws_IMAGEBOX::' 'THE FIRST INPUT TO FUNCTION gOSI_cal MUST BE A VECTOR!'],'Wrong input type');
    end