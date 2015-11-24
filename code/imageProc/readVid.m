% readvid.m: read avi file
clear all; close all; clc;

fignum=1;
dr=dir('video.avi'); % get list of movies
Nm=length(dr);       % number of movies
info=cell(1,Nm);     % set space for movie info
for n=1:Nm           % get info for each movie
  info{n}=aviinfo(dr(n).name);
end
% select a frame near the middle of movie
raw=aviread(dr(1).name,fix(info{1}.NumFrames/2));
im=double(raw.cdata(:,:,1)); % rgb all same for grayscale movie
simage(im);  % display image in false color
title('Image of an Animal to be Tracked.');
xlabel(sprintf('Figure %d.',fignum)); fignum=fignum+1;

% image correction
figure(2)
raw=aviread(dr(1).name,1);   % read first image
bg=double(raw.cdata(:,:,1)); % rgb all same for grayscale movie
simage([im bg]);  % display image in false color
title('Image of an Animal to be Tracked and Background.');
xlabel(sprintf('Figure %d.',fignum)); fignum=fignum+1;

% bright background reduction
figure(3)
ci=1-im./bg;
simage(ci);  % display image in false color
caxis(.7*[-1 1]);  % set colorscale so that zero is in the middle
colorbar;          % display color bar
title('Image of an Animal to be Tracked with Background Correction.');
xlabel(sprintf('Figure %d.',fignum)); fignum=fignum+1;

% dark background correction I
figure(4)
[nn bb]=hist(ci(:),-1:.01:1);
[mx p1]=max(nn);
nnt=nn;
nnt(bb<.3)=0;
[mx p2]=max(nnt);
dk=bb(p2);

semilogy(bb,nn,bb(p1),nn(p1),'r*',bb(p2),nn(p2),'ro','MarkerSize',10);
ylabel('Log of frequency of pixel value.');
xlabel(['Pixel value. ' sprintf('Figure %d.',fignum)]); fignum=fignum+1;
title('Histogram of Corrected Image');

% Dark background correction II
figure(5)
ci=ci/dk;
simage(ci);  % display image in false color
caxis([-.1 1.1]);  % set colorscale so that zero is blue
colorbar;          % display color bar
title('Image of an Animal to be Tracked with Background Correction.');
xlabel(sprintf('Figure %d.',fignum)); fignum=fignum+1;

% Zoom in I
figure(6)
bi=ci>.5;
simage(bi);
title('Binary Image of Animal');
xlabel(sprintf('Figure %d.',fignum)); fignum=fignum+1;

% Zoom in II
figure(7)
xprof=sum(bi,2);      % sum over horizontal direction
yprof=sum(bi,1);      % sum over vertical direction

% find position of half max of peaks (might use 0.8 of peak to get a better
% estimate of size).
xlim=find(diff(xprof>max(xprof)/2));  % x limits of animal
ylim=find(diff(yprof>max(yprof)/2));  % y limits of animal

cpx=round(mean(xlim));    % crude x position of animal
cpy=round(mean(ylim));    % crude y position of animal

sx=diff(xlim);     % crude x size (smaller than actual animal)
sy=diff(ylim);     % crude y size

plot(1:length(xprof),xprof,...
     1:length(yprof),yprof,...
     xlim,max(xprof)/2*[1 1],...
     ylim,max(yprof)/2*[1 1],...
     cpx,max(xprof)/2,'k*',...
     cpy,max(yprof)/2,'k*');

title('Identification of Center and Width of Horizontal/Vertical Profiles');
ylabel('Profile');
xlabel(['Position (pixels). ' sprintf('Figure %d.',fignum)]); fignum=fignum+1;

% Zoom in III
figure(8)
im=ci(cpx+(-sx:sx),cpy+(-sy:sy)); % examine small section
im(im<-.1)=-.1;      % remove footprints
simage(im);          % display image in false color
hold on;
th=[0:.01:2*pi 0];
D=23;[Np px py]=findcircles(double(im>.5),D,.1,1.1*D,6);
plot(py(1)+D/2*sin(th),px(1)+D/2*cos(th),'w');
plot(py(2)+D/2*sin(th),px(2)+D/2*cos(th),'w');
D=50;[Np px py]=findcircles(im,D,1,1.15*D,5);
plot(py(1)+D/2*sin(th),px(1)+D/2*cos(th),'w');
plot(py(2)+D/2*sin(th),px(2)+D/2*cos(th),'w');
hold off;
caxis([-.1 1.1]);
title('Corrected Image Zoomed in');
xlabel(sprintf('Figure %d.',fignum)); fignum=fignum+1;

% Looking for circle
figure(9)
D=50;
[x y]=ndgrid(-fix(D/2)-1:fix(D/2)+1,-fix(D/2)-1:fix(D/2)+1); % ideal particle image grid
r=abs(x+i*y);

ip=ipf(r,D,1);  % Create circle

simage([ip ip+.5*(1-ipf(r,1.15*D,1))]);
title('Test Circle');
xlabel(sprintf('Figure %d.',fignum)); fignum=fignum+1;

% Least square fit function
W=ones(size(ip));                   % Weighting factor 1.
ichi=1./chiimg(im,ip,W,[],'same');  % The inverse of the least-squares fit
                                    % function is used since it is easier
                                    % to see peaks than valleys.
[Np px py]=findpeaks(ichi,1,4,0);   % Find peaks (maxima)
[px ii]=sort(px(1:2));    % Head first (see below)
py=py(ii);

% create zoomed images
figure(10)
lx=-fix(sx/6)-1:fix(sx/6)+1;
ly=-fix(sy/6)-1:fix(sy/6)+1;
llx=(-sx:sx)/6;
lly=(-sy:sy)/6;

zi1=interp2(ly,lx,ichi(px(1)+lx,py(1)+ly),lly,llx','nearest');
zi2=interp2(ly,lx,ichi(px(2)+lx,py(2)+ly),lly,llx','nearest');

simage([ichi.*(ichi>2)+2*im.*(ichi<2) ichi; zi1 zi2]);
caxis([0 max(ichi(:))]);
hold on;
plot(py(1)+D/2*sin(th),px(1)+D/2*cos(th),'k');
plot(py(2)+D/2*sin(th),px(2)+D/2*cos(th),'k');
hold off;
title('Inverse of \chi^2 (Least-Squares Fit Function)');
colorbar;
xlabel(sprintf('Figure %d.',fignum)); fignum=fignum+1;

% Weighted fit function
figure(11)
W=ipf(r,1.15*D,1);                  % Circular weighting factor 15% larger
ichi=1./chiimg(im,ip,W,[],'same');  % The inverse of the least-squares fit
                                    % function is used since it is easier
                                    % to see peaks than valleys.
[Np px py]=findpeaks(ichi,1,4,0);   % Find peaks (maxima)
[px ii]=sort(px(1:2));    % Head first (see below)
py=py(ii);

% create zoomed images
zi1=interp2(ly,lx,ichi(px(1)+lx,py(1)+ly),lly,llx','nearest');
zi2=interp2(ly,lx,ichi(px(2)+lx,py(2)+ly),lly,llx','nearest');

simage([ichi.*(ichi>2)+2*im.*(ichi<2) ichi; zi1 zi2]);
caxis([0 max(ichi(:))]);
hold on;
plot(py(1)+D/2*sin(th),px(1)+D/2*cos(th),'k');
plot(py(2)+D/2*sin(th),px(2)+D/2*cos(th),'k');
hold off;
Title('Inverse of \chi^2 (Least-Squares Fit Function)');
colorbar;
xlabel(sprintf('Figure %d.',fignum)); fignum=fignum+1;

% use binary image
bi=ci>.5;
im=double(bi(cpx+(-sx:sx),cpy+(-sy:sy))); % examine small section

W=ipf(r,1.15*D,1);                  % Circular weighting factor 15% larger
ichi=1./chiimg(im,ip,W,[],'same');  % The inverse of the least-squares fit
                                    % function is used since it is easier
                                    % to see peaks than valleys.
[Np px py]=findpeaks(ichi,1,4,0);   % Find peaks (maxima)
[px ii]=sort(px(1:2));    % Head first (see below)
py=py(ii);

% create zoomed images
figure(12)
zi1=interp2(ly,lx,ichi(px(1)+lx,py(1)+ly),lly,llx','nearest');
zi2=interp2(ly,lx,ichi(px(2)+lx,py(2)+ly),lly,llx','nearest');

simage([ichi.*(ichi>2)+2*im.*(ichi<2) ichi; zi1 zi2]);
caxis([0 max(ichi(:))]);
hold on;
plot(py(1)+D/2*sin(th),px(1)+D/2*cos(th),'k');
plot(py(2)+D/2*sin(th),px(2)+D/2*cos(th),'k');
hold off;
Title('Inverse of \chi^2 (Least-Squares Fit Function)');
colorbar;
xlabel(sprintf('Figure %d.',fignum)); fignum=fignum+1;

% extract peak position
figure(13)
MinSep=0;        % minimum separation between peaks
Cutoff=6;        % minimum peak intensity
[Np spx spy]=findpeaks(ichi,1,Cutoff,MinSep);

simage(im);
hold on;
plot(spy,spx,'w*');
hold off;
Title('Peaks are shown in white');
xlabel(sprintf('Figure %d.',fignum)); fignum=fignum+1;