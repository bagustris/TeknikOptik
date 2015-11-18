%  result.m: show result of image tracking
clear all; close all; clc

trbifn='trackfile_bi_BN6_0_1_2_run_soft_xvid.mat';
trcifn='trackfile_ci_BN6_0_1_2_run_soft_xvid.mat';

fignum=1;

if(~exist(trbifn,'file'))  % check for trackfiles
  clear all;           % Start from scratch
  clear all
  UseBinImg=true;      % use binary image
  trackfile;
end
if(~exist(trcifn,'file')) % check for trackfiles
  clear all;           % Start from scratch
  UseBinImg=False;     % use full image
  trackfile;
end

load(trcifn);  % Load full image results
pxc=px;
pyc=py;
ppxc=ppx;
ppyc=ppy;

load(trbifn);  % Load binary image results

figure(1)
plot(px');
title('Animal Position vs. Time');
ylabel('Vertical Position (pixels)');
xlabel(['Time (frame number). ' sprintf('Figure %d.',fignum)]); fignum=fignum+1;
text(310,px(1,300),'\leftarrow Tail','HorizontalAlignment','left');
text(300,px(2,300),'Head \rightarrow','HorizontalAlignment','right');

% result 1: Sorting heads from tail
[spx ii]=sort(px); % Sort by vertical position.
spy=zeros(size(spx));
for nf=1:length(ii)    % Sort horizontal position using vertical index ii.
  spy(:,nf)=py(ii(:,nf),nf);
end

figure(2)
plot(spx');
title('Animal Position vs. Time');
ylabel('Vertical Position (pixels)');
xlabel(['Time (frame number). ' sprintf('Figure %d.',fignum)]); fignum=fignum+1;
text(310,px(1,300),'\leftarrow Tail','HorizontalAlignment','left');
text(300,px(2,300),'Head \rightarrow','HorizontalAlignment','right');

% Result 2: interpolation
[sppx ii]=sort(px+ppx); % Sort by vertical position.
                        % px (py) is list of pixel accurate positions. ppx
                        % (ppy) is a correction to px (py).  So px+ppx is
                        % sub-pixel accurate.
sppy=zeros(size(spx));
for nf=1:length(ii)    % Sort horizontal position using vertical index ii.
  sppy(:,nf)=py(ii(:,nf),nf)+ppy(ii(:,nf),nf);
end

figure(3)
plot(1:Nf,spx','.-',1:Nf,sppx');
axis([100 120 spx(1,120) spx(1,100)])
set(gca,'ygrid','on','Ytick',spx(1,120):spx(1,100))
title('Animal Position vs. Time');
ylabel('Vertical Position (pixels)');
xlabel(['Time (frame number). ' sprintf('Figure %d.',fignum)]); fignum=fignum+1;

% Result 3: Animal track
good=find(Npf==2);  % find good frames

figure(4)
plot(sppy(:,good)',sppx(:,good)');
axis('equal');
title('Animal Track');
ylabel('Vertical Position (pixels)');
xlabel(['Horizontal Position (pixels). ' sprintf('Figure %d.',fignum)]); fignum=fignum+1;

% Result 4: Animal Track
good=find(Npf==2);  % find good frames
figure(5)
plot(sppy(:,good)',sppx(:,good)','.-');
title('Animal Track');
ylabel('Vertical Position (pixels)');
xlabel(['Horizontal Position (pixels). ' sprintf('Figure %d.',fignum)]); fignum=fignum+1;

% Result 5: Binary Vs Full
[sppxc ii]=sort(pxc+ppxc); % Sort by vertical position.
                           % px (py) is list of pixel accurate positions. ppx
                           % (ppy) is a correction to px (py).  So px+ppx is
                           % sub-pixel accurate.
sppyc=zeros(size(sppxc));
for nf=1:length(ii)    % Sort horizontal position using vertical index ii.
  sppyc(:,nf)=pyc(ii(:,nf),nf)+ppyc(ii(:,nf),nf);
end
figure(6)
plot(sppy(1,good)'-10,sppx(1,good)',sppyc(1,good)',sppxc(1,good)');
title('Head Track [binary blue (offset), full green]');
ylabel('Vertical Position (pixels)');
xlabel(['Horizontal Position (pixels). ' sprintf('Figure %d.',fignum)]); fignum=fignum+1;