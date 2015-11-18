function out=zerofill(in,dNx,dNy)
% zerofill    Add zeros around an image.  
%
% Usage: out=zerofill(in,dNx,dNy)     
% Adds zeros around an image so that the final size is size(in)+[dNx dNy].

% revision history:
% 09/25/94 Mark D. Shattuck <mds> zerofill.m  

[Nx Ny]=size(in);

NNx=Nx+dNx;
NNy=Ny+dNy;
out=zeros(NNx,NNy);

out((-fix(Nx/2):ceil(Nx/2)-1)+fix(NNx/2)+1,...
    (-fix(Ny/2):ceil(Ny/2)-1)+fix(NNy/2)+1)=in;