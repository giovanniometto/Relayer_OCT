function ISOSguess = getISOSguess(im, RPE, ILM, OCTmask)

% Author: Giovanni Ometto
% Work address: C274 Tait Building City, University of London, London, EC1V 0HB (UK)
% email: giovanni.ometto@city.ac.uk
% Website: http://www.city.ac.uk
% Feb 2017; Last revision: Jun 2017

% convert image to grayscale if not already
if size(im,3)>1, greyIm = rgb2gray(im); else, greyIm = im; end

% standardisation coefficient
% coeff = size(im,1)/500;
coeff = 1;

% set some constant for the vertical-profile search
step = round(3*coeff);

% gaussian filtering
greyImG = imgaussfilt(greyIm,round(7*coeff));

% get the indexes of the top border of the mask
[~, X_nonZidx] = ind2sub(size(OCTmask),find(OCTmask));
xStart = min(X_nonZidx);
xEnd = max(X_nonZidx);

% define the steps from the first to last nonzero column of the mask
xStepsM = xStart : step : xEnd;
% if xStepsM(end) < size(im,2), xStepsM = [ xStepsM size(im,2) ]; end
xSteps = xStepsM;
ySteps = zeros(size(xStepsM));

% get the gradient based edges Y profile @ point xSteps(i) 
gradientImT = dGaussY(greyIm,round(6*coeff),round(1*coeff),1);

% get a gamma distribution function with mode @ 17.568 (calculated from
% images) starting from 10 pixels (delay) above the RPE
delay = 5; 
b=15; % scale
a=(17.568-delay)/b+1; % shape
gammaD = [ zeros(1, delay-1), gampdf((0:1:100), a, b )];

% set the gamma dist function above every reading of RPE 
idxY = flipud(repmat(1:size(im,1),size(im,2),1)');
RPEBallW = zeros(1,size(greyImG,2)); RPEBallW(RPE(1,:)) = RPE(2,:);
RPEByM = repmat(RPEBallW,size(greyImG,1),1);
idxGammaM = -round(flipud(idxY)-(RPEByM));
idxGammaM(idxGammaM>length(gammaD)) = 1;
idxGammaM(idxGammaM<1) = 1;
filtGammaM = gammaD(idxGammaM);

ILM(ILM<1) = 1; 
for i = 1:length(xSteps)
   
    % get the vertical gradient profile as the mean of column of width 2*wsize @ point X 
    halfW = round(8*coeff);
    gradYprofImT =getRowWiseColumnAvg(double(gradientImT).*filtGammaM,[xStepsM(i),1],halfW);
    
        
    % find Top RPE as first gradient peak above high intensity Y
    OptionsPks.nFirst=2; 
    OptionsPks.Treshold=.7; 
    OptionsPks.returnSortedByLocation=true;
    OptionsPks.Order=1;  
    try
        [~, L ] = getPeaks(gradYprofImT( round(ILM(2,xSteps(i))+round(15*coeff) : RPE(2,xSteps(i)) ) ), OptionsPks);
    catch
        L =[];
    end
    % manage case with no peaks
    if isempty(L), ySteps(i)=0;
    else
        ySteps(i) = L(1)+ILM(2,xSteps(i))+15*coeff;
    end


end

% % ADD first / last
nonZeroIdxs = find(ySteps);
if ySteps(1) == 0, ySteps(1) = ySteps(nonZeroIdxs(1)); xSteps(1) = 1; end 
if ySteps(end) == 0, ySteps(end) = ySteps(nonZeroIdxs(end)); xSteps(end) = size(im,2); end 

% REMOVE undetected
ySteps( xSteps == 0 ) = [];
xSteps( xSteps == 0 ) = [];

% REMOVE zeros
xSteps( ySteps == 0 ) = [];
ySteps( ySteps == 0 ) = [];

% if missing, add first and last using the gradient of the closest points
diffyT = diff(ySteps);
diffxT = diff(xSteps);
xStart = 1; xEnd = size(im,2);
if xSteps(1) ~= xStart, xSteps = [xStart xSteps]; ySteps = [ySteps(1) - diffyT(1)/diffxT(1) * (xSteps(1)-1)  ySteps] ; end
if xSteps(end) < xEnd, xSteps = [xSteps xEnd]; ySteps = [ ySteps  ySteps(end) + diffyT(1)/diffxT(1) * (size(im,2)-xSteps(end))] ; end

% interpolate results (initial guess)
ISOSguess = [interp1(xSteps,ySteps,xStart:xEnd)' (xStart:xEnd)' ];