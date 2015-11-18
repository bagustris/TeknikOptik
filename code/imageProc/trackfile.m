%% <trackfile.m> Mark D. Shattuck 7/14/2008
% Animal tracking demonstration

%% User inputs and Definitions: 
%(These may need to be changed if the animal or the image size changes.)

plotit=1;    % Set to 1 for visual feedback

D=50;      % Diameter of test circle
w=1;       % Width circles edge 
M=1.15*D;  % Diameter of weight function

MinSep=0;         % minimum separation between chi peaks 
Cutoff=6;         % minimum chi peak intensity

Profile_CutOff=60;   % minimum profile peak to track

if(~exist('UseBinImg','var')) % Can set variable outside of script.
  UseBinImg=true;      % True==use binary image--False use full image
                       % Binary is usual better if there are uncorrected
                       % spatial variations in the background.
end

res=.1;    % Interpolation resolution


%% Prepare image files:
dr=dir('video.avi'); % get list of movies
Nm=length(dr);       % number of movies
info=cell(1,Nm);     % set space for movie info
for n=1:Nm           % get info for each movie
  info{n}=aviinfo(dr(n).name);
end

Nf=info{1}.NumFrames;  % short for number of frames.
Nx=info{1}.Height;     % image size
Ny=info{1}.Width;  
ifn=dr(1).name;
[pth,base,ext,ver]=fileparts(ifn);


%% Prepare background image:
raw=aviread(ifn,1);          % read first image for bright background
bg=double(raw.cdata(:,:,1)); % rgb all same for grayscale movie

raw=aviread(ifn,fix(Nf/2));  % read frame w/ animal
im=double(raw.cdata(:,:,1)); % rgb all same for grayscale movie

ci=1-im./bg;   % Apply bright background correction.

[nn bb]=hist(ci(:),-1:.01:1);  % calculate histogram
nn(bb<.3)=0;                   % zero everything below 0.3
[mx pp]=max(nn);               % find position of peak
dk=bb(pp);                     % now any image can be normalized using:
                               % ci=(1-im./bg)/dk;
                               

%% Calculate test image
[x y]=ndgrid(-fix(D/2)-1:fix(D/2)+1,-fix(D/2)-1:fix(D/2)+1); % test image grid
r=abs(x+i*y);
 
ip=ipf(r,D,w);            % Create circle
W=ipf(r,M,1);             % Create circular weighting factor
                               

%% setup storage for results
Npf=zeros(1,Nf);  % Number of peaks per frame: should be 2 for good frames
px=zeros(2,Nf);   % Peak x positions for head and abdomen
py=zeros(2,Nf);   % Peak y positions for head and abdomen
xlim=zeros(1,2);  
ylim=zeros(1,2);  


%% Set up interpolation
[llx lly]=ndgrid(-3:res:3,-3:res:3);
ll=-3:3;
ppx=zeros(2,Nf);   % Interpolation correction
ppy=zeros(2,Nf);   % px+ppx has better resolution

%% main loop through images
for nf=1:Nf

  raw=aviread(ifn,nf);  % read frame 
  img=double(raw.cdata(:,:,1)); % rgb all same for grayscale movie
  ci=(1-img./bg)/dk;            % normalize image
  

%% Calculate crude size and position
  bi=ci>.5;                    % binary image
  xprof=sum(bi,2);      % sum over horizontal direction
  yprof=sum(bi,1);      % sum over vertical direction

  mxp=max(xprof);       % maximum of profile
  myp=max(yprof);

  mprof=max(mxp,myp);  % Max size of profile

  % find position of half max of peaks (might use 0.8 of peak to get a better
  % estimate of size).
  xlim(1)=max([find(diff(xprof>mxp/2)==1,1) 1]);        % x limits of animal
  xlim(2)=min([find(diff(xprof>mxp/2)==-1,1,'last') Nx]);
  ylim(1)=max([find(diff(yprof>myp/2)==1,1) 1]);        % y limits of animal
  ylim(2)=min([find(diff(yprof>myp/2)==-1,1,'last') Ny]);
  
  % use mprof and xylim to determine if animal is in frame
  if ((mprof > Profile_CutOff) && ...
      (xlim(1)>1) &&...
      (ylim(1)>1) &&...
      (xlim(2)<Nx) &&...
      (ylim(2)<Ny))  % else no need to track
    % Finish calculating crude size and position
  
    cpx=round(mean(xlim));    % crude x position of animal
    cpy=round(mean(ylim));    % crude y position of animal

    sx1=diff(xlim);     % crude x size (smaller than actual animal)
    sy1=diff(ylim);     % crude y size
 

%% deal with zooming in at boundaries
    sx2=sx1;
    sy2=sy1;
    if (sx1>cpx-1); sx1=cpx-1; end   % don't go over lower bound
    if (sy1>cpy-1); sy1=cpy-1; end
    if (sx2>Nx-cpx); sx2=Nx-cpx; end   % don't go over upper bound
    if (sy2>Ny-cpy); sy2=Ny-cpy; end   

    if (UseBinImg)
      im=double(bi(cpx+(-sx1:sx2),cpy+(-sy1:sy2))); % examine small section
    else
      im=ci(cpx+(-sx1:sx2),cpy+(-sy1:sy2)); % examine small section
      im(im<-.1)=-.1;                       % remove footprints
    end

%% begin tracking frame nf
    ichi=1./chiimg(im,ip,W,[],'same');  % The inverse of the least-squares fit

    [Np spx spy]=findpeaks(ichi,1,Cutoff,MinSep);  % find peaks
    Npf(nf)=min(2,Np);        % Check for at least 2 peaks
    % Store first 2 peak positions
    if (Npf(nf)>0)
      % findpeaks returns peak sorted by peak height so first peaks are the
      % largest peaks.  Also convert to full image coordinates 
      px(:,nf)=spx(1:Npf(nf))+cpx-sx1-1;  
      py(:,nf)=spy(1:Npf(nf))+cpy-sy1-1;  
      for np=1:Npf(nf);
        iichi=interp2(ll,ll,ichi(spx(np)+ll,spy(np)+ll),lly,llx,'spline');
        [mx pp]=max(iichi(:));
        ppx(np,nf)=llx(pp);
        ppy(np,nf)=lly(pp);
      end
    end

    
%% Visual feedback
    if (plotit==1)
      simage(im);
      hold on;
      plot(py(:,nf)-cpy+sy1+1,px(:,nf)-cpx+sx1+1,'w*');
      hold off;
      drawnow;
    end
  end
  fprintf('Processing Frame %04d/%04d. mprof=%d. xlim=[%d,%d]. ylim=[%d,%d]. %d points found.\n',nf,Nf,mprof,xlim,ylim,Npf(nf));
end


%% Save Data
if(UseBinImg)  % Tag file name with type of image analyzed.
  IType='bi';
else
  IType='ci';
end
save(['trackfile_' IType '_' base]);

