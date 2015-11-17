% intDif.m : demo of interference and difraction
clear
clc
thetamax=pi/50;
while(1)
    a=input('Enter slit width (in micro meter): ');
    a=a*1e-6;
    d=input('Enter slit seperation (in mm): ');
    d=d*1e-3;
    l=input('Enter wavelength (in nm): ');
    l=l*1e-9;
    s=input('Enter slit to screen distance (in m): ');
    theta=-thetamax:1e-5:thetamax;
    y=s*tan(theta);
    alpha=pi*a*sin(theta)/l;
    beta=pi*d*sin(theta)/l;
    x1=cos(beta).^2;            % Interference term
    x2=(sin(alpha)./alpha).^2;  % Diffraction term
    x=x1.*x2;                   % Combined effect
    plot(y,x,'b',y,x2,'--r');
    title('Double slit diffraction');
    xlabel('Distance in m');
    ylabel('Intensity');
    hold all;
    ch= input('Press 1 to continue and 0 to exit: ');
if ch == 0
    break;
end
end