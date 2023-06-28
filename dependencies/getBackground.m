function bkgrdROI=getBackground(avg,bw0,xy,dilateN,width,height)
[m,n,ks]=size(bw0);
bwAll_outer=false(m,n);
bwAll_cell=false(m,n);
I=cell(ks,1);
bwROIAll=false(m,n,ks);
bwAll_back=false(m,n);
for ii=1:ks
    bw=bw0(:,:,ii);
    
    bwAll_cell=or(bw,bwAll_cell);
    bw2=bwmorph(bw,'dilate',dilateN);
    bwAll_outer=or(bw2,bwAll_outer);
end
for ii=1:ks
    bw=bw0(:,:,ii);
    cx=xy{ii}(:,1);
    cy=xy{ii}(:,2);
    t=regionprops(bw,'centroid');
    x0=t.Centroid(1);
    y0=t.Centroid(2);
    [x1,x2,y1,y2]=getSquare(x0,y0,avg,width,height);
    bwROI=false(m,n);
    bwROI(y1:y2,x1:x2)=true;
    bwROI(bwAll_outer)=false;
    bwAll_back=or(bwAll_back,bwROI);
    I{ii}=find(bwROI(:));
    bwROIAll(:,:,ii)=bwROI;
end
bkgrdROI.I=I;
bkgrdROI.width=width;
bkgrdROI.height=height;
bkgrdROI.bw=bw0;
bkgrdROI.avg=avg;
bkgrdROI.xy=xy;
bkgrdROI.bwAll_back=bwAll_back;
bkgrdROI.bwAll_cell=bwAll_cell;
bkgrdROI.bwAll_outer=bwAll_outer;
bkgrdROI.bwROIAll=bwROIAll;
function [x1,x2,y1,y2]=getSquare(x0,y0,img,width,height)
marginLeft=0;
marginTop=0;
marginRight=0;
marginBottom=0;
x1=x0-width/2;
y1=y0-height/2;
x1=max([marginLeft+1,x1]);
y1=max([marginTop+1,y1]);

x2=x1+width-1;
y2=y1+height-1;
[m,n]=size(img);
if x2>n-marginRight
    x2=n-marginRight;
end
if y2>m-marginBottom
    y2=m-marginBottom;
end
x1=round(x1);
x2=round(x2);
y1=round(y1);
y2=round(y2);