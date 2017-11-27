function yEdgesIm = dGaussY (im, width, sigma, whiteToBlack)

% FILTEDGEY = returns the image filtered with a filter of size 2*width to
% enhance edges in the Y dimension. If whiteToBlack is true will enhance 
% edges of transitions from bright to dark pixels along increasing y  
% indexes (i. e. downwards) - dark to white otherwise, which is the default 
% option.
%
% Author: Giovanni Ometto
% Work address: C274 Tait Building City, University of London, London, EC1V 0HB (UK)
% email: giovanni.ometto@city.ac.uk
% Website: http://www.city.ac.uk
% Feb 2017; Last revision:

% convert image to grayscale if not already
if size(im,3)>1, greyIm = rgb2gray(im); else, greyIm = im; end

if nargin<4
    whiteToBlack = 0;
end

% get the gaussian filter
t = (-width:width);
ssq = sigma^2;
gau = exp(-(t.*t)/(2*ssq))/(2*pi*ssq);

% get the 2d gaussian-gradient filter 
[x,y]=meshgrid(-width:width,-width:width);
filtEdgeY=(2*-whiteToBlack+1)*y.*exp(-(y.*y+x.*x)/(2*ssq))/(pi*ssq);

% smooth gray image and filter across columns
imSmoothX=imfilter(greyIm,gau,'conv','replicate');   % run the filter across rows
imSmooth=imfilter(imSmoothX,gau','conv','replicate'); % smooth image across columns
yEdgesIm = imfilter(imSmooth, filtEdgeY, 'conv','replicate');