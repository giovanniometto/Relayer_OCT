function ISOS = snakeISOS(im, ISOSguess)

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
diffyL = ISOSguess(2*cw,1)-ISOSguess(1,1);
diffyR = ISOSguess(end-2*cw,1)-ISOSguess(end,1);
diffx = 2*cw;

ISOSguessPadL = flipud( repmat(ISOSguess(1,1),cw,1) - diffyL/diffx * (1:cw)' );
ISOSguessPadR = repmat(ISOSguess(end,1),cw,1) - diffyR/diffx * (1:cw)' ;

ISOSguess = [ISOSguessPadL; ISOSguess(:,1); ISOSguessPadR] ;
ISOSguess = [ISOSguess (1:length(ISOSguess))'];

% standardisation coefficient
% coeff = size(im,1)/500;
coeff = 1;

% Set parameters and calculate the snake
Options.Wline=0; 
Options.Wedge=10; 
Options.Sigma1=round(.5*coeff); 
Options.Sigma2=round(.8*coeff);
Options.Alpha=1; %100?
Options.Beta=10; % 10?
Options.Gamma=1;
Options.Kappa=.4;
Options.Iterations=50;
% Options.Verbose=true;

ISOS = Snake1D(greyIm,ISOSguess,Options);

% exclude the padding
ISOS = ISOS(:,cw+1:end-cw);
ISOS(1,:) = 1:size(ISOS,2);