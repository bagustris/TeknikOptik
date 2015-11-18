function [Np px py ichi]=findcircles(im,D,w,M,CutOff)
% findcircles Find best matching circular objects in image
% Usage: [Np px py ichi]=findcircles(im,D,w,M,CutOff)
%
% Returns Np centers at px, py of the best matches (chi>CutOff) to a circle
% with diameter D and edge thickness w.  Any pixels outside of a circle of
% diameter M are ignored  Can return 1/chi map as ichi. 

% revision history:
% 07/18/08 Mark D. Shattuck <mds> findcircles.m
%          combination of chiimg and findpeaks.

[x y]=ndgrid(-fix(M/2)-1:fix(M/2)+1,-fix(M/2)-1:fix(M/2)+1); % ideal particle image grid
r=abs(x+i*y);
ichi=(1./chiimg(im,ipf(r,D,w),ipf(r,M,1),[],'same'));
[Np px py]=findpeaks(ichi,1,CutOff,0);