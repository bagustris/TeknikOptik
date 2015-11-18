function h=simage(img,numcolor,axistype,minimg,maximg);
% simage Display scaled image.
% Usage: h=simage(img,numcolor,axistype,minimg,maximg);
%
% Displays img scaled to use the full colormap.  numcolor is only included
% for backwards compatibility [use colormap].  axistype is passed to axis.
% The default is 'image'.  min[max]img can be used to manually set the max
% and min values of the scaled image. 

% revision history:
% 05/12/94 Mark D. Shattuck <mds> simage.m  
% 07/22/95 mds add min[max]img
% 04/20/00 mds changed to use imagesc [numcolor now irrelevant]

imscale=0;
if ~exist('axistype')
  axistype=['image'];
end

if ~exist('minimg') & ~exist('maximg')
  imscale=1;
else
  if (~exist('maximg'))  
    maximg=max(img(:));
  end
  if (~exist('minimg'))  
    minimg=min(img(:));
  end
  if (maximg==0 & minimg==0)
    minimg=min(img(:));
    maximg=max(img(:));
  end
end

if ~exist('numcolor')
  numcolor=64;
end


if imscale==0
  h=imagesc(img,[minimg maximg]);
else
 h=imagesc(img);
end

axis(axistype);

