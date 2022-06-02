function [data3D,varargout]=groupZprojected_registered(data,avgNumber)
    if avgNumber==1
        data3D=data;
    else
        [m,n,p]=size(data);
        k1=floor(p/avgNumber);
        k2=mod(p,avgNumber);
        data3D=zeros(m,n,k1);
        if k2==0
            i1=cell(k1,1);
            for ii=1:k1
                i1{ii}=[(ii-1)*avgNumber+1:(ii-1)*avgNumber+avgNumber];
            end
        else
            i1=cell(k1+1,1);
            for ii=1:k1
                i1{ii}=[(ii-1)*avgNumber+1:(ii-1)*avgNumber+avgNumber];
            end   
            i1{k1+1}=[(k1+1-1)*avgNumber+1:p];
        end
        k1=length(i1);
        options.maxLoop=15; % maximum iterations allowed to register images; 
        options.minSSE=2; % threshold of SSE; if SSE is smaller than this value, then ierations stop;    
        options.previousRefFlag=0; % when it is 0, register each image with resp
        for ii=1:k1
            x=data(:,:,i1{ii});
            x=imgRegistration(x,options);
            data3D(:,:,ii)=mean(x,3);
            data(:,:,i1{ii})=x;

        end        
    end

  
    options.previousRefFlag=1;
    [data3D,yxShiftAll]=imgRegistration(data3D,options);
    if nargout>1 && avgNumber>1
        for ii=1:k1
            x=data(:,:,i1{ii});
            yxShift=yxShiftAll(ii,:);
            for jj=1:size(x,3);
                x(:,:,jj)=circshift(x(:,:,jj),yxShift);
            end
            data(:,:,i1{ii})=x;
        end
        varargout{1}=data;
    end