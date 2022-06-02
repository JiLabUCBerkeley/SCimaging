function data=readSingleTif(file,method,input)
if nargin==1
    method=1;
    input.evenShift=0;
elseif nargin==2
    input.evenShift=0;
end

if method==1
    info=imfinfo(file);
p=length(info);
n=info(1).Width;
m=info(1).Height;

data=zeros(m,n,p,'uint16');
for ii=1:p
    img=(imread(file, ii, 'Info', info));
    img(2:2:end,:)=circshift(img(2:2:end,:),[0,input.evenShift]);
    data(:,:,ii)=img;
end
elseif method==2
    InfoImage=imfinfo(file);
    mImage=InfoImage(1).Width;
    nImage=InfoImage(1).Height;
    NumberImages=length(InfoImage);
    data=zeros(nImage,mImage,NumberImages,'uint16');
    FileID = tifflib('open',file,'r');
    rps = tifflib('getField',FileID,Tiff.TagID.RowsPerStrip);

    for i=1:NumberImages
       tifflib('setDirectory',FileID,i);
       % Go through each strip of data.
       rps = min(rps,nImage);
       for r = 1:rps:nImage
          row_inds = r:min(nImage,r+rps-1);
          stripNum = tifflib('computeStrip',FileID,r);
          img = tifflib('readEncodedStrip',FileID,stripNum);
          img(2:2:end,:)=circshift(img(2:2:end,:),[0,input.evenShift]);
          data(row_inds,:,i)=img;
       end
    end
    tifflib('close',FileID);
end