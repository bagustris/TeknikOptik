function [cr,over]=pgrid(spx,spy,NNx,NNy,rect,Np,D,rad)  
% pgrid  Create a local grid <cr=cx+i*cy> centered on each particle and an overlap matrix <over> 
% Usage: [cr,over]=pgrid(spx,spy,NNx,NNy,rect,Np,D,rad)  
%
% Creates a local grid [see ndgrid] for the Voronoi volume of each Np point
% centered at the position specified by spx and spy.  The union of the
% grids should cover and NNx x NNy image inside the rectangle [rect].  The
% maximum size of an individual grid is 2*ceil(D/2).  If rad=1 then return
% abs(cr) instead.  over is an image whose value is the index of spx/spy
% [1:Np] at each pixel in the Voronoi volume of each point. 

% revision history:
% 01/16/01 Mark D. Shattuck <mds> calcimg.m  
% 02/30/06 mds rename pgrid.m return cr and over instead of ci
% 03/26/14 mds (green)


if (isempty(spx));
  cr=zeros(NNx,NNy);
  over=cr;
  return;
end
yl=rect(1);
yh=rect(2);
xl=rect(3);
xh=rect(4);
[xxx yyy]=ndgrid(1:NNx,1:NNy);
cr=max(NNx,NNy)*ones(NNx,NNy);
over=zeros(NNx,NNy);
lll=-ceil(D/2):ceil(D/2);
kkk=lll;
for np=1:Np
  ll=round(spx(np)) + lll;
  kk=round(spy(np)) + kkk;
  ll=ll((ll<=xh)&(ll>=xl));
  kk=kk((kk<=yh)&(kk>=yl));
  if(numel(ll) && numel(kk))       % particles influence must be inside rect
    if(rad)
      X=abs(xxx(ll,kk)-spx(np)+1i*(yyy(ll,kk)-spy(np)));
      [nx,~]=size(X);
      [ix iy]=find(cr(ll,kk)>=X);
      over(ix+ll(1)-1+NNx*(iy+kk(1)-2))=np;
      cr(ix+ll(1)-1+NNx*(iy+kk(1)-2))=X(ix+nx*(iy-1));
    else
      X=xxx(ll,kk)-spx(np)+1i*(yyy(ll,kk)-spy(np));
      [nx,~]=size(X);
      [ix iy]=find(abs(cr(ll,kk))>=abs(X));
      over(ix+ll(1)-1+NNx*(iy+kk(1)-2))=np;
      cr(ix+ll(1)-1+NNx*(iy+kk(1)-2))=X(ix+nx*(iy-1));
    end
  end
end