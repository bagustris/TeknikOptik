function [chiimg Wip2]=chiimg(img,ip,W,Wip2,range)
% chiimg     Calculate chi-squared image
% Usage: [chiimg Wip2]=chiimg(img,ip,W,Wip2,range);
%
% Calculates an image of chi-squared of the form 
% chiimg=int(W(x-x0)(img(x)-ip(x-x0))^2 dx) using convolution.  Chi-squared
% is an image which can be larger than (range=='full'[default]), smaller
% than (range=='valid'), or the same size as (range='same') the input image
% img. chiimg is minimum where img and the test image ip are most alike in
% a squared-difference sense.  W ([default]W==ip) limits the area to be
% consider by weighting chiimg.  

% revision history:
% 08/04/00 Mark D. Shattuck <mds> chiimg.m  
% 01/30/04 mds added return Wip2
% 02/22/04 mds added range option
% 09/22/07 mds update for non symmetric W and Ip

if(~exist('range','var'))
  range='full';
end
if(~exist('W','var') || isempty(W))
  W=ip;
end
if(~exist('Wip2','var') || isempty(Wip2))  % Wip2 can be pre calculated since it does not depend on img
  blk=ones(size(img));                     % Blank image
  Wip2=(conv2(blk,ip.^2.*W,range));        % Weighting factor
end

ip=ip(end:-1:1,end:-1:1);  % Flip for convolution
W=W(end:-1:1,end:-1:1);    % Flip for convolution

chiimg=1+(-2*conv2(img,ip.*W,range)+conv2(img.^2,W,range))./Wip2;    % best fit ignoring overlap  
