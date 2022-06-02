function  data=addOrientation_v2(data,angleLeft2,angleNO,framePerSti,options)
% v2, make it compact having an options
if nargin==4
        options.preMoving=2;
        options.postMoving=0;
        options.blank=0;
end
 p=size(data,3);% p is the frame number in total
 %% creating gratings  
    gratingSize=20;% size of grating
    gratingSize_padding=4*gratingSize; % bigger template so that when we rotate the grating, there is no black;
%     beta=0;
%     fshift=4/n2;
    movingSpeed=1;
%     beta=beta-pi/2;
%     mapTmp2=zeros(gratingSize,gratingSize);
%     signFlag=1;
%     for ii=1:gratingSize
%         if mod(ii,1/(fshift*2))==0
%             signFlag=signFlag*-1;
%         
%         end
%         if signFlag==1
%             
%         else
%             mapTmp2(ii,:)=1;
%         end
%         
%     end
    gratingSize_x=1:gratingSize_padding;gratingSize_x=gratingSize_x(:);
    fshift=3/gratingSize;% period per pixel; 3 periods in total;
    cosFun=cos(2*pi*fshift*gratingSize_x)*ones(1,gratingSize_padding);
    mapTmp2=cosFun>=0;

    figure(1);clf;imshow(mapTmp2,[])
    x1=gratingSize_padding/2-gratingSize/2;x1=floor(x1);
    x2=x1+gratingSize-1;
    y1=x1;y2=x2;
    

    t=1;
    ii2=0;
    mapTmp2=uint16(mapTmp2);
for ss=1:floor(p/angleNO/framePerSti)% floor(p/angleNO/framePerSti) is repetition number
    for ii=1:angleNO
        dx=0;
        for jj=1:framePerSti
            img2=data(:,:,t);
             beta=angleLeft2(t);
                if (jj>=1 && jj<= options.preMoving) % static period before grating moves
                    dx=0;
                    dx=floor(dx);
                    mapTmp3=circshift(mapTmp2,[-dx,0]);
                    mapTmp4=imrotate(mapTmp3,beta,'nearest');
                    mapTmp5=mapTmp4(y1:y2,x1:x2);

                    if jj==1 && ii2==0
                        imgMax=max(img2(:))*0.7;
                        ii2=1;
                    end
                    mapTmp5=mapTmp5*(imgMax);

                        if options.blank==1 % if the value is 1, make the image blank
                            mapTmp5=imgMax*0.5*(0*mapTmp5+1);
                        end
                        

                    img2(1:gratingSize,end-gratingSize+1:end)=mapTmp5;     
                elseif jj>options.preMoving && jj<=framePerSti-options.postMoving % moving grating
                    dx=dx+movingSpeed;
                    dx=floor(dx);
                    mapTmp3=circshift(mapTmp2,[-dx,0]);
                    mapTmp4=imrotate(mapTmp3,beta,'nearest');
                    mapTmp5=mapTmp4(y1:y2,x1:x2);
                    if jj==1 && ii2==0
                        imgMax=mean(max(img2))*0.7;
                        ii2=1;
                    end
                    mapTmp5=mapTmp5*imgMax;

                    img2(1:gratingSize,end-gratingSize+1:end)=mapTmp5;   
                elseif jj>framePerSti-options.postMoving% static period after grating moves
                        if options.blank==1 % if the value is 1, make the image blank
                            mapTmp5=imgMax*0.5*(0*mapTmp5+1);
                        end                    
                    img2(1:gratingSize,end-gratingSize+1:end)=mapTmp5;   
                    
                    %%
           
                end
           
            data(:,:,t)=img2;
             t=t+1;
        end
    end
end
trash=0;