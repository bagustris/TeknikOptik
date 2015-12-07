% Script to calculate field distribution from a point source% ---------------------------------------% Specify the dimensions of the system% ---------------------------------------% Default values: Wavelength = 4.0e-3 m  % Distance source to screen = 5.0 cm = 5.0e-2 m% Range of detection on the screen = +/- 1.2 cm = 2.4e-2 m  lambda   = 4.0e-3;  scrnDist = 5.0e-2;  scrnWdth = 2.4e-2;% (xs,ys) is position of the source% In this exercise the source is at the origin  xs = 0;    ys = 0;  % The amplitude of the field at the source is arbitrary   A  = 1;   % ------------------------------------------% Specify the detector points% ------------------------------------------% Let detector points be a rectangular array of points % between the source and the screen.% Divide range of x and y into N segments (lengths dX and dY)  N       = 5;                            % default N=5, change to 500  dX      = scrnDist/N;  Xcoords = [dX : dX : scrnDist]';  dY      = scrnWdth/N;  Ycoords = [-scrnWdth/2 : dY : +scrnWdth/2]';% Coordinates of the detector points form an N+1 by N matrix  [xd,yd] = meshgrid(Xcoords,Ycoords); % ------------------------------------------% Calculate fields at detector points at t=0% ------------------------------------------% Distances between source and detectors  r = sqrt((xd-xs).^2 + (yd-ys).^2); % Value of the field at the points (xd,yd)  E0 = A * cos(2*pi*r/lambda) ./r;% plotfigure(1); mesh(xd,yd,E0)figure(2); ContourPlot(xd,yd,E0);