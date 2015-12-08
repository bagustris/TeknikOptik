% <trackframe.m> Mark D. Shattuck 3/29/2008
% Particle tracking demonstration

% User Inputs
D=12;          % Initial Diameter Guess
w=1.3;         % Initial Width Guess
Cutoff=5;      % minimum peak intensity
MinSep=5;      % minimum separation between peaks

hi=250;  % hi and lo values come the image histogram
lo=10;   % hi/lo=typical pixel value outside/inside 

maxDwnr=10;       % maximum number of calls to cidDw Newton solver.  
mindelD=.0001;    % minimum change in D before stopping

maxnr=5;       % maximum number of calls to cidp2 Newton solver.  
mindelchi2=1;  % minimum change in chi2 before stopping

% setup for ideal particle

ss=2*fix(D/2+4*w/2)-1;         % size of ideal particle image
os=(ss-1)/2;                   % (size-1)/2 of ideal particle image
[xx yy]=ndgrid(-os:os,-os:os);  % ideal particle image grid
r=abs(xx+i*yy);    % radial coordinate

raw=imread('test.bmp');  % load image
[Nx Ny]=size(raw);       % image size

ri=(hi-double(raw))/(hi-lo);  % normalize image 

% find pixel accurate centers using chi-squared

[Np px py]=findpeaks(1./chiimg(ri,ipf(r,D,w)),1,Cutoff,MinSep);  % find maxima

% Minimizing chi-squared for sub-pixel accuracy

[cxy over]=pgrid(px-os,py-os,Nx,Ny,[1 Nx 1 Ny],Np,2*os+3,0); % create local grid centered on each particle and overlap matrix
ci=ipf(cxy,D,w);  % create calculated image

di=ci-ri;                             % Calculate difference image 
chi2=sum(di(:).^2);                   % Calculate Chi-Squared
fprintf('Chi-Squared=%f\n',chi2);

chi2o=chi2;di0=di;  % save for later

% find best D and w
nr=0;
delD=1e99;
while((abs(delD)>mindelD) && (nr<maxDwnr))
  [delD delw]=cidDw(abs(cxy),di,D,w);
  D=D+delD;
  w=w+delw;
  ci=ipf(abs(cxy),D,w);
  di=(ci-ri);
  fprintf('.');
  nr=nr+1;
end
fprintf('\n');
chi2=sum(di(:).^2);
fprintf('Chi-Squared=%f\n',chi2);

% Find best positions
nr=0;
delchi2=1e99;
while((abs(delchi2)>mindelchi2) && (nr<maxnr))
  [dpx dpy]=cidp2(cxy,over,di,Np,D,w);  %
  px=px+dpx;
  py=py+dpy;
  [cxy over]=pgrid(px-os,py-os,Nx,Ny,[1 Nx 1 Ny],Np,2*os+3,0); % create local grid centered on each particle and overlap matrix
  ci=ipf(cxy,D,w);  % create calculated image
  di=(ci-ri);
  delchi2=chi2-sum(di(:).^2); % Calculate change in Chi-Squared
  chi2=chi2-delchi2;          % Calculate Chi-Squared
  fprintf('.');
  nr=nr+1;
end
fprintf('\n');
chi2=sum(di(:).^2);                     % Calculate Chi-Squared
fprintf('Chi-Squared=%f\n',chi2);

h=figure(2); set(h,'Position',[100 100 600 400],'Color',[1 1 1]);
simage([di.^2 di0.^2]); caxis([0 .25]); % Show the chi-squared images
title(sprintf('New \\chi^2=%6.2f           Original \\chi^2=%6.2f',chi2,chi2o),'fontsize',15)
