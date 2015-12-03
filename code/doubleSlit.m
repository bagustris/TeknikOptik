% intDif.m : demo of interference and difraction

clear
clc
thetamax=pi/50;

a=input('Enter slit width (in micro meter): '); % ex=30
a=a*1e-6;
d=input('Enter slit seperation (in mm): ');     % ex=0.15
d=d*1e-3;
l=input('Enter wavelength (in nm): ');          % ex=557
l=l*1e-9;
s=input('Enter slit to screen distance (in m): ');  % ex= 0.5
    
theta=-thetamax:1e-5:thetamax;
y=s*tan(theta);
alpha=pi*a*sin(theta)/l;
beta=pi*d*sin(theta)/l;
x1=cos(beta).^2;            % Interference term
x2=(sin(alpha)./alpha).^2;  % Diffraction term
x=x1.*x2;                   % Combined effect
    
plot(x,y,'b',x2,y,'--r');
title('Double slit diffraction');
xlabel('Intensity');
ylabel('Distance in m');
hold all;