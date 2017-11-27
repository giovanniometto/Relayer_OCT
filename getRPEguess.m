function RPEguess  = getRPEguess(im, ILM, OCTmask)

% Author: Giovanni Ometto
% Work address: C274 Tait Building City, University of London, London, EC1V 0HB (UK)
% email: giovanni.ometto@city.ac.uk
% Website: http://www.city.ac.uk
% Jan 2017; Last revision: Jun 2017

if nargin<4
    flagplot =0;
end

% standardisation coefficient
% coeff = size(im,1)/500;
coeff = 1;

% convert image to grayscale if not already
if size(im,3)>1, greyIm = rgb2gray(im); else, greyIm = im; end

% set some constant for the vertical-profile search
step = round(4*coeff);
halfW = round(8*coeff);

% gaussian filtering
greyImG = imgaussfilt(greyIm,round(9*coeff));

% get the gradient based edges Y profile @ point xSteps(i) 
gradientImB = dGaussY(greyIm,round(7*coeff),round(.7*coeff),0);
% figure(4); imshow(gradientImB,[]);

% get the indexes of the top border of the mask
[~, X_nonZidx] = ind2sub(size(OCTmask),find(OCTmask));
xStart = min(X_nonZidx);
xEnd = max(X_nonZidx);

% define the steps from the first to last nonzero column of the mask
xStepsM = xStart : step : xEnd;
% if xStepsM(end) < size(im,2), xStepsM = [ xStepsM size(im,2) ]; end
xSteps = xStepsM;
ySteps = zeros(size(xStepsM));

% ******************************** RPE ************************************

% Detection of the centre Y of hyper-reflective complex thorugh vertical
% inspection of the intensity profile
for i = 1:length(xSteps)
    
    % get the vertical intensity profile as the mean of column of width 2*halfW @ point xStepsM 
    smoothYprofIm = getRowWiseColumnAvg(greyImG,[xSteps(i),1],halfW);
    OptionsPks.nFirst=0; 
    OptionsPks.Treshold=.5; 
    OptionsPks.returnSortedByLocation=false;
    OptionsPks.Order=-1;  
    [~, L ] = getPeaks(smoothYprofIm, OptionsPks);
    L(L < ILM(2,xSteps(i))+20*coeff) = [];
    if isempty(L), yStepsM(i)=0;  
    else
        yStepsM(i) = L(1);
    end
    
    % get the vertical gradient profile as the mean of column of width 2*halfW @ point xStepsM 
    halfW = round(10*coeff);
    gradYprofImB = getRowWiseColumnAvg(gradientImB,[xStepsM(i),1],halfW); 
%     gradYprofImB = getRowWiseColumnAvg(gradientImB,[xStepsM(i),1],halfW); 

    if yStepsM(i)~=0 &&  yStepsM(i) > ILM(2,xSteps(i))+20*coeff+3
    
        % find Bottom RPE as first gradient paek below high intensity Y
        OptionsPks.nFirst=1; 
        OptionsPks.Treshold=0.5; 
        OptionsPks.returnSortedByLocation=true;
        OptionsPks.Order=-1;  
        [~, L ] = getPeaks(round(gradYprofImB( max(1,yStepsM(i)-round(1.3*coeff)) : min(size(gradYprofImB,1),yStepsM(i)+round(13*coeff)) )), OptionsPks);
        % manage case with no peaks
        if isempty(L), ySteps(i)=0;
        else
            ySteps(i) = L+yStepsM(i); 
        end
        
    else
        ySteps(i)=0;
    end

end

% % ADD first / last
nonZeroIdxs = find(ySteps);
if ySteps(1) == 0, ySteps(1) = ySteps(nonZeroIdxs(1)); xSteps(1) = 1; end 
if ySteps(end) == 0, ySteps(end) = ySteps(nonZeroIdxs(end)); xSteps(end) = size(im,2); end 
% if yStepsB(1) == 0
%     xStepsB(1) = 1;
%     yStepsB = yStepsB(nonZeroIdxs(1)) - (yStepsB(nonZeroIdxs(2))-yStepsB(nonZeroIdxs(1)))/(xStepsB(nonZeroIdxs(2))-xStepsB(nonZeroIdxs(1))) * (xStepsB(nonZeroIdxs(1))-1) ; 
% end 
% if yStepsB(end) == 0
%     xStepsB(end) = size(im,2);
%     yStepsB = yStepsB(nonZeroIdxs(end)) - (yStepsB(nonZeroIdxs(end))+yStepsB(nonZeroIdxs(end-1)))/(xStepsB(nonZeroIdxs(end))-xStepsB(nonZeroIdxs(end-1))) * (xStepsB(nonZeroIdxs(end))-xStepsB(nonZeroIdxs(end-1))) ; 
% end 


% REMOVE undetected
ySteps( xSteps == 0 ) = [];
xSteps( xSteps == 0 ) = [];

% REMOVE zeros
xSteps( ySteps == 0 ) = [];
ySteps( ySteps == 0 ) = [];

% % REMOVE high gradients measurements
heightCoeff = (500/size(im,1));
gradB = abs(diff(ySteps)./diff(xSteps)).*heightCoeff;
while sum(gradB > 1) ~= 0
    ySteps( gradB > 1) = [];
    xSteps( gradB > 1) = [];
    gradB = abs(diff(ySteps)./diff(xSteps))*heightCoeff;
end

% if missing, add first and last using the gradient of the closest points
diffyB = diff(ySteps);
diffxB = diff(xSteps);
xStart = 1; xEnd = size(im,2);
if xSteps(1) ~= xStart,  ySteps = [ySteps(1) - diffyB(1)/diffxB(1) * (xSteps(1)-1)  ySteps] ;  xSteps = [xStart xSteps]; end
if xSteps(end) < xEnd, xSteps = [xSteps xEnd]; ySteps = [ ySteps  ySteps(end) + diffyB(1)/diffxB(1) * (size(im,2)-xSteps(end))] ; end

% interpolate results (initial guess)
RPEguess = [interp1(xSteps,ySteps,xStart:xEnd)' (xStart:xEnd)' ];
