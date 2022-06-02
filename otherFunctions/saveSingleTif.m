function saveSingleTif(file, data,method)

if nargin==2
    method=1;
end
[m,n,p]=size(data);
if m~=n
    method=2;
elseif n<=64
    method=2;
end

if method==1
    maxP=4000;
    if p>maxP
        p2=floor(p/maxP);
        indP=cell(p2,1);
        for ii=1:p2-1
            indP{ii}=(ii-1)*maxP+(1:maxP);
        end
        indP{end}=((p2-1)*maxP+1):p;
        for ii=1:p2
            file2=strrep(file,'.tif',[num2str(ii,'%03d'),'.tif']);
            saveSingleTif_method1(file2,data(:,:,indP{ii}));
        end
    else
        
        saveSingleTif_method1(file,data);
    end
       
elseif method==2
    maxP=4000;
    if p>maxP
        p2=floor(p/maxP);
        indP=cell(p2,1);
        for ii=1:p2-1
            indP{ii}=(ii-1)*maxP+(1:maxP);
        end
        indP{end}=((p2-1)*maxP+1):p;
        for ii=1:p2
            file2=strrep(file,'.tif',[num2str(ii,'%03d'),'.tif']);
            saveSingleTif_method2(file2,data(:,:,indP{ii}));
        end
    else
        
        saveSingleTif_method2(file,data);
    end
    
end
function saveSingleTif_method1(file,data)
if exist(file)==0
else
    delete(file)
end  
    tifstreamobj2 = TifStream(file, size(data,2), size(data,1), 16);
    for jj=1:size(data,3)
            tifstreamobj2.appendFrame(data(:,:,jj));
    end 
    tifstreamobj2.close; 
function saveSingleTif_method2(file,data)
if exist(file)==0
else
    delete(file)
end    
data=uint16(data);
    bigtiff=false;
    fname=file;
    bitspersamp=16;
     if bigtiff
        t = Tiff(fname,'w8');
    else
        t = Tiff(fname,'w');
    end
    tagstruct.ImageLength = size(data,1);
    tagstruct.ImageWidth = size(data,2);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    if bitspersamp==16
        tagstruct.BitsPerSample = 16;
    end
    if bitspersamp==32
        tagstruct.BitsPerSample = 32;
    end
    tagstruct.SamplesPerPixel = 1;
    tagstruct.RowsPerStrip = 256;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Software = 'MATLAB';
    t.setTag(tagstruct);
    t.write(data(:,:,1));
    numframes = size(data,3);
    divider = 10^(floor(log10(numframes))-1);
    tic
    for i=2:numframes
        t.writeDirectory();
        t.setTag(tagstruct);
        t.write(data(:,:,i));
        if (round(i/divider)==i/divider)
            fprintf('Frame %d written in %.0f seconds, %2d percent complete, time left=%.0f seconds \n', ...
                i, toc, i/numframes*100, (numframes - i)/(i/toc));
        end
    end
    t.close();    

    