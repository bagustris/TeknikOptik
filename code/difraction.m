% difraction.m
% The theory of optics predicts that the diffraction pattern produced
% by a plane wave incident on an optical mask with a small aperture
% is described, at a distance, by the Fourier transform of the mask.

% create logical array
clear all; close all; clc;
n = 2^10;
M = zeros(n);				% mask

I = 1:n;
x = I-n/2;
y = n/2-I;
[X,Y] = meshgrid(x,y);
R = 10;                     % radius
A = (X.^2 + Y.^2 <= R^2);   % circular apperture
M(A) = 1;

figure(1)
imagesc(M)
colormap([0 0 0; 1 1 1])
axis image
title('{\bf Circular Aperture}')

% Use fft2 to compute the two-dimensional DFT of the mask
% fftshift to rearrange the output 
% so that the zero-frequency component is at the center.


D1 = fft2(M);
D2 = fftshift(D1);

figure(2)
imagesc(abs(D2))
axis image
colormap(hot)
title('{\bf Diffraction Pattern}')

% Use to bring out details of the DFT in regions 
% where the amplitude is small

D3 = log2(D2);

figure(3)
imagesc(abs(D3))
axis image
colormap(hot)
title('{\bf Enhanced Diffraction Pattern}')
