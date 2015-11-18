function [dpx,dpy,dipxx]=cidp2(cxy,over,di,Np,D,w)  
% cidp2    Calculate one Newton's step toward minimizing di^2 over particle centers. 
% Usage: [dpx,dpy]=cidp2(cxy,over,di,Np,D,w)  
%
% Calculates change in px and py needed to move di^2 closer to a minimum.  

% revision history:
% 09/14/00 Mark D. Shattuck <mds> cidp2.m  
% 02/15/02 mds implement over and separate calc per particle
% 06/25/04 mds limit change to maxdr=2;
% 04/30/07 mds changed meaning of w to 1/w
% 03/26/14 mds made -sign consistant with cidDw.m (green)

w=1/w;      % w is real width
maxdr=20;    % max change

A=zeros(2,2);
dp=zeros(Np,2);

[~, idx]=sort(over(:));  % create list of pixels for each particle
nn=histc(over(:),0:Np);
cnn=cumsum(nn);  

rr=abs(cxy)+eps;         % useful numbers
rr3=rr.*rr.*rr+eps;
xx=real(cxy);
yy=imag(cxy);
xx2=xx.^2;
yy2=yy.^2;
tanh1=tanh((rr-D/2)*w);
sech2=sech((rr-D/2)*w).^2;
 
dipx=w.*xx.*sech2./2./rr;  % first derivatives
dipy=w.*yy.*sech2./2./rr;
dipxx=w.*sech2.*(2*w*xx2.*rr.*tanh1-yy2)./2./rr3;  % second derivatives
dipyy=w.*sech2.*(2*w*yy2.*rr.*tanh1-xx2)./2./rr3;
dipxy=w.*xx.*yy.*sech2.*(2*w*rr.*tanh1+1)./2./rr3;

chix=di.*dipx;    % first derivatives
chiy=di.*dipy;
chixx=dipx.^2+di.*dipxx;   % second derivatives
chiyy=dipy.^2+di.*dipyy;
chixy=dipx.*dipy+di.*dipxy;

for np=1:Np  % loop over particles
  ii=idx(cnn(np)+1:cnn(np+1));
  b=[sum(chix(ii)) sum(chiy(ii))];
  A(1,1)=sum(chixx(ii));
  A(1,2)=sum(chixy(ii));
  A(2,1)=A(1,2);
  A(2,2)=sum(chiyy(ii));
  dp(np,:)=-b/A;   % Newton's step
end

dpx=dp(:,1);
dpy=dp(:,2);
dpc=dpx+1i*dpy;
ii=find(abs(dpc)>maxdr); % max a unit step only for large steps
dpx(ii)=dpx(ii)./abs(dpc(ii));
dpy(ii)=dpy(ii)./abs(dpc(ii));

