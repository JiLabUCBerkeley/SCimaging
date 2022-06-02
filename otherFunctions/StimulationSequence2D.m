function [order,point_Temp_Spatial_trial]=StimulationSequence2D(filePath,fileName,pp)

    
    
p=length(fileName);
order=cell(p,1);
point_Temp_Spatial_trial=zeros(p,4);
id=1:11;
if length(pp)==1
    pp=ones(p,1)*pp;
end
for ii=1:p
    file1=fileName(ii).name;
    fileTmp=dir(fullfile(filePath,[file1(id),'*.mat']));
    if ~isempty(fileTmp)
    [order{ii},point_Temp_Spatial_trial(ii,:)]=StimulationSequence_singleMat(fullfile(filePath,fileTmp(1).name),pp(ii));
 
    end

   
    
end
function [order,point_Temp_Spatial_trial]=StimulationSequence_singleMat(file,pp)
    point_Temp_Spatial__trial=zeros(1,4);
    matObj=matfile(file);
if ~isempty(whos(matObj,'contrastRaw2'))

    
    trialNO=matObj.trialNO;
    freq_pool=matObj.freq_pool;
    contrastRaw2=matObj.contrastRaw2;
    angleRaw2=matObj.angleRaw2;
    pSection=length(freq_pool);
    contrastAll=zeros(pSection,1);
    angleAll=zeros(pSection,1);
    m1=matObj.m1;n1=matObj.n1;
    for ii=1:pSection
        [II,JJ]=ind2sub([m1,n1],freq_pool(ii));
        contrast=contrastRaw2(II,JJ);
        contrastAll(ii)=contrast;
        angleAll(ii)=angleRaw2(II,JJ);
    end
    contrastRaw=matObj.contrastRaw;
    angleRaw=matObj.angleRaw;
    spatialP=length(contrastRaw);
    temporalP=length(angleRaw);
    
    
    point_Temp_Spatial_trial=[pp/temporalP/spatialP/trialNO,temporalP,spatialP,trialNO];
    

    point=point_Temp_Spatial_trial(1);
    temporalP=point_Temp_Spatial_trial(2);
    spatialP=point_Temp_Spatial_trial(3);
    trialNO=point_Temp_Spatial_trial(4);
    order.p1=point;
    order.p2=temporalP;
    order.p3=spatialP;
    order.p4=trialNO;
    order.pAll=[order.p1,order.p2,order.p3,order.p4];
    order1=reshape((1:point).'*ones(1,temporalP*spatialP*trialNO),pp,1);    
    order.order1=order1;    
    order2b=ones(point,1)*angleAll.';
    order2=order2b(:);
    order3b=ones(point,1)*contrastAll.';
    order3=order3b(:);    
    order.order4=reshape(ones(point*temporalP*spatialP,1)*(1:trialNO),pp,1);
    order.order1=order1;
    order.order2=order2;
    order.order3=order3;
    order.orderAll=[order.order1,order.order2,order.order3,order.order4];
    order.pAll=[order.p1,order.p2,order.p3,order.p4];


else


    spatial_freq_actual=matObj.spatial_freq_actual;
    freq_pool=matObj.freq_pool;
    freq_pool=freq_pool/max(freq_pool)*max(matObj.spatial_freq_actual);
    cyclespersec_pool=matObj.cyclespersec_pool;
    trialNO=matObj.trialNO;
    temporalP=length(cyclespersec_pool);
    spatialP=length(spatial_freq_actual);

    point_Temp_Spatial_trial=[pp/temporalP/spatialP/trialNO,temporalP,spatialP,trialNO];
    point=point_Temp_Spatial_trial(1);
    temporalP=point_Temp_Spatial_trial(2);
    spatialP=point_Temp_Spatial_trial(3);
    trialNO=point_Temp_Spatial_trial(4);
    temporalOrder1=ones(point,1)*cyclespersec_pool.';
    % temporalOrder1=temporalOrder1(:);
    temporalOrder2=temporalOrder1(:)*ones(1,spatialP*trialNO);
    temporalOrder=temporalOrder2(:);
    order.temporalOrder=temporalOrder;

    spatialOrder1=ones(point*temporalP,1)*(freq_pool(:)).';
    spatialOrder=spatialOrder1(:);
    order.spatialOrder=spatialOrder;
    order.p1=point;
    order.p2=temporalP;
    order.p3=spatialP;
    order.p4=trialNO;
    order.pAll=[order.p1,order.p2,order.p3,order.p4];

    order1=reshape((1:point).'*ones(1,temporalP*spatialP*trialNO),pp,1);
    order.order1=order1;
    order.order2=temporalOrder;
    order.order3=spatialOrder;
    order.order4=reshape(ones(point*temporalP*spatialP,1)*(1:trialNO),pp,1);
    order.orderAll=[order.order1,order.order2,order.order3,order.order4];

end
% point_Temp_Spatial__trial=zeros(1,4);
% matObj=matfile(file);
% 
% spatial_freq_actual=matObj.spatial_freq_actual;
% freq_pool=matObj.freq_pool;
% freq_pool=freq_pool/max(freq_pool)*max(matObj.spatial_freq_actual);
% cyclespersec_pool=matObj.cyclespersec_pool;
% trialNO=matObj.trialNO;
% temporalP=length(cyclespersec_pool);
% spatialP=length(spatial_freq_actual);
% 
% point_Temp_Spatial_trial=[pp/temporalP/spatialP/trialNO,temporalP,spatialP,trialNO];
% point=point_Temp_Spatial_trial(1);
% temporalP=point_Temp_Spatial_trial(2);
% spatialP=point_Temp_Spatial_trial(3);
% trialNO=point_Temp_Spatial_trial(4);
% temporalOrder1=ones(point,1)*cyclespersec_pool.';
% % temporalOrder1=temporalOrder1(:);
% temporalOrder2=temporalOrder1(:)*ones(1,spatialP*trialNO);
% temporalOrder=temporalOrder2(:);
% order.temporalOrder=temporalOrder;
% 
% spatialOrder1=ones(point*temporalP,1)*(freq_pool(:)).';
% spatialOrder=spatialOrder1(:);
% order.spatialOrder=spatialOrder;
% order.p1=point;
% order.p2=temporalP;
% order.p3=spatialP;
% order.p4=trialNO;
% order.pAll=[order.p1,order.p2,order.p3,order.p4];
% 
% order1=reshape((1:point).'*ones(1,temporalP*spatialP*trialNO),pp,1);
% order.order1=order1;
% order.order2=temporalOrder;
% order.order3=spatialOrder;
% order.order4=reshape(ones(point*temporalP*spatialP,1)*(1:trialNO),pp,1);
% order.orderAll=[order.order1,order.order2,order.order3,order.order4];









