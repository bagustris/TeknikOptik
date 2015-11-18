function [Npf,spx,spy]=findpeaks(img,mk,CutOff,MinSep)
% findpeaks     Find intensity peaks in a image
% Usage: [Npf,spx,spy]=findpeaks(img,mask,CutOff,MinSep)
%
% Finds all pixels in img*mask, which are larger than their 8 nearest
% neighbors and have intensities greater than 'Cutoff' and are separated
% from all other peaks by at least 'MinSep' pixels.  
%
% revision history:
% 02/24/01 Mark D. Shattuck <mds> chiimg.m  
% 03/21/03 mds added mask
% 04/03/03 mds change implementation of mask
% 10/02/07 mds added MinSep
% 07/18/08 mds added sort by peak height;


[NNx NNy]=size(img);

uu=mk;
for n=-1:1
  for m=-1:1
    if(~(m==0 && n==0))  
      uu=(img>img(rem((1:NNx)+NNx+n-1,NNx)+1,rem((1:NNy)+NNy+m-1,NNy)+1)) & uu;
    end
  end
end

g1=find((img(:).*uu(:))>CutOff);
Npf=length(g1);
[spx spy]=ind2sub([NNx NNy],g1);
[junk iii]=sort(img(g1));
Xn=repmat(spx(iii)+i*spy(iii),1,Npf);
dd=abs(Xn.'-Xn);
iix=iii((sum(tril(dd<MinSep & dd~=0)))~=0);
g1=g1(setdiff(1:Npf,iix));
[junk iii]=sort(img(g1),'descend');
Npf=length(g1);
[spx spy]=ind2sub([NNx NNy],g1(iii));


