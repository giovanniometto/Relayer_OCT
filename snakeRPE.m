function [RPE] = snakeRPE(im,RPEguess)

% Author: Giovanni Ometto
% Work address: C274 Tait Building City, University of London, London, EC1V 0HB (UK)
% email: giovanni.ometto@city.ac.uk
% Website: http://www.city.ac.uk
% Feb 2017; Last revision: Jun 2017

% convert image to grayscale if not already
if size(im,3)>1, greyIm = rgb2gray(im); else, greyIm = im; end

% add size(im,2)/20 black padding to the left and right of the image
h = size(greyIm,1);
cw = round(size(greyIm,2)/20);
blackColumn = zeros(h, cw);
greyIm = [ blackColumn greyIm blackColumn ];

% extend guess to the width of the padded image
% if missing, add first and last using the gradient of the closest points
diffyL = RPEguess(2*cw,1)-RPEguess(1,1);
diffyR = RPEguess(end-2*cw,1)-RPEguess(end,1);
diffx = 2*cw;

RPEguessPadL = flipud( repmat(RPEguess(1,1),cw,1) - diffyL/diffx * (1:cw)') ;
RPEguessPadR = repmat(RPEguess(end,1),cw,1) - diffyR/diffx * (1:cw)' ;

% % Just repeat the last value
% RPEguessPadL = repmat(RPEguess(1,1),cw,1);
% RPEguessPadR = repmat(RPEguess(end,1),cw,1);

RPEguess = [RPEguessPadL; RPEguess(:,1); RPEguessPadR] ;
RPEguess = [RPEguess (1:length(RPEguess))'];

% standardisation coefficient
% coeff = h/500;
coeff = 1; 

% Set parameters and calculate the snake
Options.Wline=0; 
Options.Gamma=1;
% Options.Sigma2=round(2*coeff);
% Options.Sigma1=round(.6*coeff); 
Options.Sigma1=round(2*coeff); 
Options.Sigma2=round(3*coeff);
Options.Wedge=10;

Options.Alpha= 1;   % 0;
Options.Beta= 10^4; % 10^5;
Options.Kappa=.4;

% Options.Wedge=50;
% Options.Sigma1=round(1*coeff);
% Options.Alpha=15; 
% Options.Beta=20; 
% Options.Kappa=.5;

Options.Iterations=20;
% Options.Verbose=true;

RPE = Snake1D(greyIm,RPEguess,Options);

% exclude the padding
RPE = RPE(:,cw+1:end-cw);
RPE(1,:) = 1:size(RPE,2);