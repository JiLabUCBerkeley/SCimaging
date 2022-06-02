 function  [data,yxShiftAll2]=imgRegistration_largeDataSet(data,options)
 [m, n, o]=size(data);
 avgN=40;
 p=floor(o/avgN);
 indN=zeros(o,1)+p;
 indNTmp=ones(avgN,1)*(1:p);
 indN(1:length(indNTmp(:)))=indNTmp(:);
 [data3,data3Std,angleLeft2]=avgRepeatNO(data,indN);
%  [data3,yxShift,yxShiftAll]=dtfreg_fast_2mean(data3);
%   [data3,yxShift,yxShiftAll]=dtfreg_fast_2First(data3);
  [data3,yxShiftAll]=imgRegistration(data3,options);
%  [reg_frames,yxShift,yxShiftAll]=dtfreg_fast_2First(frames2reg)
%  [yxShiftAll]= dtfreg_fast_shifts(fft2(data3),fft2(data3(:,:,1)),'corrected');
yShift=ones(avgN,1)*yxShiftAll(:,1).';
xShift=ones(avgN,1)*yxShiftAll(:,2).';
yxShiftAll2=[yShift(:),xShift(:)];
%  yxShiftAll2=removeMean(yxShiftAll2);
o1=size(yxShiftAll2,1);
 for ii=1:o1
     if isnan(yxShiftAll2(ii,:))
     else
         data(:,:,ii)=circshift(data(:,:,ii),yxShiftAll2(ii,:));
     end
     
 end
 if o1<o
     for ii=o1+1:o
         if isnan(yxShiftAll2(end,:))
         else
             data(:,:,ii)=circshift(data(:,:,ii),yxShiftAll2(end,:));
         end         
     end
 end
     function shift=removeMean(shift)
         [m,n]=size(shift);
         for ii=1:n
             [N,X]=hist(shift(:,ii),20);
             [Y,I]=max(N);
             shift(:,ii)=shift(:,ii)-round(X(I));
         end