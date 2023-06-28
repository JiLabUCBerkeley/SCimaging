function [gDSI,pref_direction,exception]=gDSI_cal_updated(ydata,datapoints,options)
% A direction-selectivity index (DSI) as the vector sum of responses
% normalized by the scalar sum of responses (such that the index varies
% between 0 and 1).
% The angle of this vector sum defined the preferred direction of each
% cell.
%created by Wenzhi Sun, July 27 2015

if nargin==2
    options.negativeMethod=2;
end

if isvector(ydata)
    if eq(numel(ydata),datapoints)
        if iscolumn(ydata), ydata=ydata'; end
        theta=0:(2*pi/datapoints):(2*pi-0.1);
        %for 8 datapoints-measurement theta=[0 pi/4 pi/2 pi3/4 pi pi5/4 pi3/2 and pi7/4]
        if options.negativeMethod==2
            if min(ydata)<0, ydata=ydata-min(ydata); end %remove the negtive value
        elseif options.negativeMethod==1
            ydata(ydata<0)=0;
        end
        

        gDSI=abs(sum(ydata.*exp(1i.*theta))/sum(ydata));
        pref_direction=angle(sum(ydata.*exp(1i.*theta)));
        %P = angle(Z) returns the phase angles, in radians, for each element of complex array Z. The angles lie between ±pi.
        %to wrap the phase angle to [0-2pi)
        if pref_direction<0, pref_direction=pref_direction+pi*2; end
        exception=[];
    else
        gDSI=nan;
        pref_direction=nan;
        exception=MException(['ws_IMAGEBOX::' 'THE DATA POINTS OF INPUT DATA DOSN''T MATCH THE INPUT DATAPOINTS!'],'Iputs Not Match');
    end
else
    gDSI=nan;
    pref_direction=nan;
    exception=MException(['ws_IMAGEBOX::' 'THE FIRST INPUT TO FUNCTION gDSI_cal MUST BE A VECTOR!'],'Wrong input type');   
end  