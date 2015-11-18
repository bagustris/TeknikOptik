function [delD delw]=cidDw(cr,di,D,w)  
% cidDw    Calculate one Newton's step toward minimizing di^2 over D and w. 
% Usage: [dpx,dpy]=cidDw(cxy,di,D,w)  
%
% Calculates change in D and w needed to move di^2 closer to a minimum.  

% revision history:
% 09/14/00 Mark D. Shattuck <mds> cidDw.m  
% 04/30/07 mds changed meaning of w to 1/w
% 03/26/14 mds fix 1/w mistake did not matter much since w is near 1

w=1/w;

A=zeros(2,2);

rp=cr-D/2;
tanh1=tanh(rp*w);
sech2=sech(rp*w).^2;

dipD=w.*sech2/4;
dipw=-rp/2.*sech2;
dipDD=w^2/4*tanh1.*sech2;
dipww=rp.^2.*tanh1.*sech2;
dipDw=sech2.*(1-2*w*rp.*tanh1)/4;

chiD=di.*dipD;
chiw=di.*dipw;
chiDD=dipD.^2+di.*dipDD;
chiww=dipw.^2+di.*dipww;
chiDw=dipD.*dipw+di.*dipDw;

b=[sum(chiD(:)) sum(chiw(:))];
A(1,1)=sum(chiDD(:));
A(1,2)=sum(chiDw(:));
A(2,1)=A(1,2);
A(2,2)=sum(chiww(:));
delDw=-b/A;
delD=delDw(1);
delw=-delDw(2)/(w+delDw(2))/w;
