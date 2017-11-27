function [ILMguess] = getILMguess(im,OCTmask)

%GETILMSNAKE1D - returns the coordinates of the ILM found by a 1D active
%contour allowing movement of the points in the vertical direction.
%
% Author: Giovanni Ometto
% Work address: C274 Tait Building City, University of London, London, EC1V 0HB (UK)
% email: giovanni.ometto@city.ac.uk
% Website: http://www.city.ac.uk
% Jan 2017; Last revision: Jun 2017

% standardisation coefficient
% coeff = size(im,1)/495;
coeff = 1; 

% convert image to grayscale if not already
if size(im,3)>1, greyIm = rgb2gray(im); else, greyIm = im; end
cols = size(im,2);
rows = size(im,1);


% create a new image and relative mask with zero padding at the top
emptyRows = round(rows*.1);
% emptyRows = round(0);
zerosM = zeros(rows+emptyRows,cols);
zerosM(emptyRows+1:end,:) = greyIm;
greyIm = zerosM;
falseM = false(rows+emptyRows,cols);
falseM(emptyRows+1:end,:) = OCTmask;
OCTmask = falseM;

% gaussian filtering
greyImG = imgaussfilt(greyIm,round(4*coeff));
% greyImG2 = imgaussfilt(greyIm,round(10*coeff));
% Initial guess of the border points of inner limiting membrane as upper 
% peak of the sobel filtered image @ points xSteps

% set some constants for the vertical-profile search
step = round(15*coeff);
halfW = round(7*coeff);

% get the gradient based edges Y profile @ point xSteps(i) 
gradientIm = imgradient(imgaussfilt(greyImG,1),'sobel'); 

% get the indexes of the top border of the mask
[~, topY] = max( OCTmask, [], 1 );
% [~, X_nonZidx] = ind2sub(size(OCTmask),find(OCTmask));
% xStart = min(X_nonZidx);
% xEnd = max(X_nonZidx);

% define the steps from the first to last nonzero column of the mask
xSteps = 1 : step : size(im,2);
if xSteps(end) ~= size(im,2), xSteps = [ xSteps size(im,2) ] ; end
ySteps = zeros(size(xSteps));

for i = 1:length(xSteps) % column-wise peak search
    
    gradYprofIm = getRowWiseColumnAvg(gradientIm,[xSteps(i),1],halfW);
    
    % smooth vertical profile to remove noise
    g = gausswin(round(step/2)); 
    smoothGradYprofIm = conv(gradYprofIm, g, 'same');
    
    % find peaks > 40% of the max peak and get the uppest Y coordinate
    OptionsPks.nFirst=4; 
    OptionsPks.Treshold=0.2; 
    OptionsPks.returnSortedByLocation=true;
    OptionsPks.Order=1;  
    [P, L ] = getPeaks(smoothGradYprofIm, OptionsPks);
    
    % manage case with no peaks
    if isempty(L)
        ySteps(i)=0; 
    elseif L(1) < topY(i)+10*coeff && P(1) < .7*max(P) % if the peak in gradient is found at the top border and is not particularly high compared to the max, take the next peak
         ySteps(i)=(L(2)); 
    else
        ySteps(i)=(L(1));  
    end
    
end

% % ADD first / last
nonZeroIdxs = find(ySteps);
if ySteps(1) == 0, ySteps(1) = ySteps(nonZeroIdxs(1)); end 
if ySteps(end) == 0, ySteps(end) = ySteps(nonZeroIdxs(end)); end 

% REMOVE undetected
ySteps( xSteps == 0 ) = [];
xSteps( xSteps == 0 ) = [];

% REMOVE readings at the top edge
[xSteps, ySteps] = removeEdgeReadings( xSteps, ySteps, topY, OCTmask, greyIm, coeff);

% REMOVE high gradient readings
[xSteps, ySteps] = removeHighGradReadings( xSteps, ySteps, OCTmask, greyIm, 1);

% interpolate results (initial guess)
ILMguess = [interp1(xSteps,ySteps,1:size(im,2))'-emptyRows (1:size(im,2))' ];

% *************************************************************************




function  [xStepsNoEdge, yStepsNoEdge] = removeEdgeReadings( xSteps, ySteps, topY, OCTmask, greyIm, coeff)

% REMOVE spikes to the top edge of the scan if the characteristic of the
% pixels below resemble those of *likely* ares of the vitreous

xStepsNoEdge = xSteps;
yStepsNoEdge = ySteps;
readingsAtTop = (yStepsNoEdge<topY(xStepsNoEdge)+4*coeff);

if sum(readingsAtTop)~=0
    
    falseReadings = false(size(readingsAtTop));
    [~,Ys] = meshgrid(1:size(greyIm,2),1:size(greyIm,1));
    
    % get a mask of pixels above the line connecting all readings but those at the top
    xStepsC = xStepsNoEdge;
    yStepsC = yStepsNoEdge;
    yStepsC(end) = yStepsNoEdge(find(~readingsAtTop,1,'last'));
    yStepsC(1) = yStepsNoEdge(find(~readingsAtTop,1,'first'));
    % remove all readings at top with the exception of first and last
    xStepsC([false readingsAtTop(2:end-1) false]) = [];
    yStepsC([false readingsAtTop(2:end-1) false]) = [];
    clearedLine = interp1(xStepsC,yStepsC,1:size(greyIm,2));
    Yc = repmat(clearedLine,size(Ys,1),1);
    maskCleared = and(Ys<Yc,OCTmask);
    
    % get mask from actual readings (likely only vitreous)
    safeVGelB = interp1(xStepsNoEdge,yStepsNoEdge,1:size(greyIm,2));
    Yvg = repmat(safeVGelB,size(Ys,1),1);
    maskVGel = and(Ys<Yvg,OCTmask);
    meanVGel = mean(greyIm(maskVGel));
    stdVGel = std(greyIm(maskVGel));

    % get the areas under the readings to the top
    maskAUTopReadings = and(maskCleared,~maskVGel);
    
    % start indexes
    idxReadingsAtTop = find(readingsAtTop);
    startIdxs = idxReadingsAtTop([true diff(find(readingsAtTop))>1])-1;
%     if startIdxs(1) == 0, startIdxs(1) = 1; end
    
    % end indexes
    idxReadingsAtTop = find(readingsAtTop);
    endIdxs = idxReadingsAtTop([ fliplr(abs(diff(fliplr(idxReadingsAtTop))))>1 true])+1;
%     if endIdxs(1) == length(endIdxs)+1, endIdxs(end) = length(endIdxs); end
    
    for i=1:length(startIdxs)
        imSector = greyIm(:, xStepsNoEdge(max(1,startIdxs(i))): xStepsNoEdge(min(endIdxs(i),length(xStepsNoEdge)) ));
        maskSector = maskAUTopReadings(:, xStepsNoEdge(max(1,startIdxs(i))): xStepsNoEdge(min(endIdxs(i),length(xStepsNoEdge)) ));
        if sum(maskSector(:))~=0
            falseReadings(startIdxs(i)+1:endIdxs(i)-1) = mean(imSector(maskSector)) < (meanVGel+stdVGel);
        end
    end
    
% REMOVE false detections at the top
yStepsNoEdge( falseReadings ) = [];
xStepsNoEdge( falseReadings ) = [];

% add last if missing
    if xStepsNoEdge(end) ~= xSteps(end)
        gradLast = (yStepsNoEdge(end)-yStepsNoEdge(end-1))/(xStepsNoEdge(end)-xStepsNoEdge(end-1));
        yEstimate =   yStepsNoEdge(end) - (xSteps(end) - xStepsNoEdge(end)) * gradLast;
        yStepsNoEdge = [ yStepsNoEdge yEstimate ];
        xStepsNoEdge = [ xStepsNoEdge xSteps(end) ];    
    end
    % add first if missing
    if xStepsNoEdge(1) ~= xSteps(1)
        gradFirst = (yStepsNoEdge(2)-yStepsNoEdge(1))/(xStepsNoEdge(2)-xStepsNoEdge(1));
        yEstimate =   yStepsNoEdge(1) - (xStepsNoEdge(1) - xSteps(1)) * gradFirst;
        yStepsNoEdge = [ yEstimate yStepsNoEdge ];
        xStepsNoEdge = [ xSteps(1) xStepsNoEdge ];
    end
end


function  [ xSteps, ySteps] = removeHighGradReadings( xSteps, ySteps, OCTmask, greyIm, treshGrad)

% REMOVE readings between high/low gradients if the characteristic of the
% pixels below resemble those of *likely* ares of the vitreous or if less
% than 3 segments and the longest is much longer than the other(s)

grad = diff(ySteps)./diff(xSteps);
highGradNeg = [(grad)>treshGrad false]; 
highGradPos = [(grad)<-treshGrad false];

if any ( or(highGradNeg,highGradPos) )
    
    xStepsNoHighGrad = xSteps; 
    yStepsNoHighGrad = ySteps;
    
    [startConnecetdIdx, endConnecetdIdx] = getSegmentsStartEnd(xSteps,ySteps,highGradPos,highGradNeg);
    boundHighGradSeg = [startConnecetdIdx endConnecetdIdx];
    segmentsN = size(boundHighGradSeg,1);
    segmentsLengths = endConnecetdIdx - startConnecetdIdx;
    [longestSegLength,idxLongestSeg] = max(segmentsLengths);
    shortestSegLengths = segmentsLengths(segmentsLengths~=longestSegLength);
%     
    % if less than 3 segments and the max is longer than five time the
    % other(s), take only the longest
    if segmentsN == 2 && (longestSegLength > 5*shortestSegLengths)
       
    RETURNxSteps = xSteps(startConnecetdIdx(idxLongestSeg):endConnecetdIdx(idxLongestSeg));
    RETURNySteps = ySteps(startConnecetdIdx(idxLongestSeg):endConnecetdIdx(idxLongestSeg));
        
    elseif segmentsN == 3 && (longestSegLength > 5*min(shortestSegLengths))
        
%         RETURNxSteps1 = xSteps(startConnecetdIdx(idxLongestSeg):endConnecetdIdx(idxLongestSeg));
%         RETURNySteps1 = ySteps(startConnecetdIdx(idxLongestSeg):endConnecetdIdx(idxLongestSeg));
        RETURNxSteps = xSteps;
        RETURNySteps = ySteps;
        if (startConnecetdIdx(idxLongestSeg) == 1)

            RETURNxSteps(startConnecetdIdx(idxLongestSeg+1):endConnecetdIdx(idxLongestSeg+1) ) = [];
            RETURNySteps(startConnecetdIdx(idxLongestSeg+1):endConnecetdIdx(idxLongestSeg+1) ) = [];
        elseif(endConnecetdIdx(idxLongestSeg) == length(xSteps)) 
            RETURNxSteps = xSteps;
            RETURNySteps = ySteps;
            RETURNxSteps(startConnecetdIdx(idxLongestSeg-1):endConnecetdIdx(idxLongestSeg-1) ) = [];
            RETURNySteps(startConnecetdIdx(idxLongestSeg-1):endConnecetdIdx(idxLongestSeg-1) ) = [];
        end
    
    else

        
        if min(find(highGradPos))>min(find(highGradNeg))    % if first is a negative high grad, highGrad segments would be the ones at even rows - remove the odd indexes to preserve them from deleting the correspondi values   
            boundHighGradSeg(2:2:segmentsN,:)=[] ;
        elseif and(startConnecetdIdx(1) == endConnecetdIdx(1), highGradNeg(1) )   % if the first point only is high grad get the even rows - remove the odd
            boundHighGradSeg(2:2:segmentsN,:)=[] ;
        else                                                % if first is a positive high grad, highGrad segments would be the ones at odd rows - remove the even   
            boundHighGradSeg(1:2:segmentsN,:)=[] ;
        end

        % get the line connecting non high-grad measurements

        for i=fliplr(1:size(boundHighGradSeg,1))
            xStepsNoHighGrad(boundHighGradSeg(i,1):boundHighGradSeg(i,2)) = [];
            yStepsNoHighGrad(boundHighGradSeg(i,1):boundHighGradSeg(i,2)) = [];
        end

        % add last if missing
        if xStepsNoHighGrad(end) ~= xSteps(end)
            gradLast = (yStepsNoHighGrad(end)-yStepsNoHighGrad(end-1))/(xStepsNoHighGrad(end)-xStepsNoHighGrad(end-1));
            yEstimate =   yStepsNoHighGrad(end) - (xSteps(end) - xStepsNoHighGrad(end)) * gradLast;
            yStepsNoHighGrad = [ yStepsNoHighGrad yEstimate ];
            xStepsNoHighGrad = [ xStepsNoHighGrad xSteps(end) ];    
        end
        % add first if missing
        if xStepsNoHighGrad(1) ~= xSteps(1)
            gradFirst = (yStepsNoHighGrad(2)-yStepsNoHighGrad(1))/(xStepsNoHighGrad(2)-xStepsNoHighGrad(1));
            yEstimate =   yStepsNoHighGrad(1) - (xStepsNoHighGrad(1) - xSteps(1)) * gradFirst;
            yStepsNoHighGrad = [ yEstimate yStepsNoHighGrad ];
            xStepsNoHighGrad = [ xSteps(1) xStepsNoHighGrad ];
        end

        falseReadings = false(size(xSteps));

        [~,Ys] = meshgrid(1:size(greyIm,2),1:size(greyIm,1));

        % get a mask of pixels above the line connecting all readings but those
        % between high gradients
        clearedLine = interp1(xStepsNoHighGrad,yStepsNoHighGrad,1:size(greyIm,2));
        Yc = repmat(clearedLine,size(Ys,1),1);
        maskCleared = and(Ys<Yc,OCTmask);

        % get mask from actual readings (likely only vitreous)
        safeVGelB = interp1(xSteps,ySteps,1:size(greyIm,2));
        Yvg = repmat(safeVGelB,size(Ys,1),1);
        maskVGel = and(Ys<Yvg,OCTmask);
        meanVGel = mean(greyIm(maskVGel));
        stdVGel = std(greyIm(maskVGel));

        % get the areas under the readings to the top
        maskAUTopReadings = and(maskCleared,~maskVGel);

        % start indexes
        startConnecetdIdx = boundHighGradSeg(:,1)-1;
    %     startConnecetdIdx(startConnecetdIdx<1) = 1;
        startIdxs = startConnecetdIdx;

        % end indexes
        endConnecetdIdx = boundHighGradSeg(:,2)+1;
    %     endConnecetdIdx(endConnecetdIdx>length(xSteps)) = length(xSteps);
        endIdxs = endConnecetdIdx;

        for i=1:length(startIdxs)
            imSector = greyIm(:, xSteps(max(1,startIdxs(i))): xSteps(min(endIdxs(i),length(xSteps)) ));
            maskSector = maskAUTopReadings(:, xSteps(max(1,startIdxs(i))): xSteps(min(endIdxs(i),length(xSteps)) ));
            if sum(maskSector(:))~=0
                falseReadings(startIdxs(i)+1:endIdxs(i)-1) = mean(imSector(maskSector)) < (meanVGel+stdVGel);
            end
        end

        % REMOVE false detections
        RETURNySteps = ySteps;
        RETURNySteps( falseReadings ) = [];
        RETURNxSteps = xSteps;
        RETURNxSteps( falseReadings ) = [];
    end
        
    % add last if missing
    if RETURNxSteps(end) ~= xSteps(end)
        gradLast = (RETURNySteps(end)-RETURNySteps(end-1))/(RETURNxSteps(end)-RETURNxSteps(end-1));
        yEstimate =   RETURNySteps(end) - (xSteps(end) - RETURNxSteps(end)) * gradLast;
        RETURNySteps = [ RETURNySteps yEstimate ];
        RETURNxSteps = [ RETURNxSteps xSteps(end) ];    
    end
    % add first if missing
    if RETURNxSteps(1) ~= xSteps(1)
        gradFirst = (RETURNySteps(2)-RETURNySteps(1))/(RETURNxSteps(2)-RETURNxSteps(1));
        yEstimate =   RETURNySteps(1) - (RETURNxSteps(1) - xSteps(1)) * gradFirst;
        RETURNySteps = [ yEstimate RETURNySteps ];
        RETURNxSteps = [ xSteps(1) RETURNxSteps ];
    end
    
    xSteps = RETURNxSteps;
    ySteps = RETURNySteps;

end


function [startConnecetdIdx, endConnecetdIdx] = getSegmentsStartEnd(xS,~,highGradPos,highGradNeg)

idxPosGrad = find(highGradPos);
idxNegGrad = find(highGradNeg); 

startConnecetdIdx = 1; 
endConnecetdIdx = [];

if isempty(idxNegGrad)
    currentHighGrad = idxPosGrad(1);
    checkPosList = false;         % check tag for the other list
elseif isempty(idxPosGrad)
    currentHighGrad = idxNegGrad(1);
    checkPosList = true;         % check tag for the other list
else
    currentHighGrad = min(idxPosGrad(1),idxNegGrad(1));     % get index of first measured high gradient
    checkPosList = min(idxNegGrad)<min(idxPosGrad);         % check tag for the other list
end

endConnecetdIdx = [ endConnecetdIdx ;...                % update segment end boundaries
          currentHighGrad ];
startConnecetdIdx = [ startConnecetdIdx; ...            % update segment start boundaries
            currentHighGrad+1 ];


if checkPosList                                         % delete first selected index from the list
    idxNegGrad(idxNegGrad==currentHighGrad) = [];
else
    idxPosGrad(idxPosGrad==currentHighGrad) = [];
end

while ~isempty(idxPosGrad) || ~isempty(idxNegGrad)
    if checkPosList  && ~isempty(idxPosGrad)                                   
        currentHighGrad = idxPosGrad(1);
        endConnecetdIdx = [ endConnecetdIdx ;...            % update segment end boundaries
              currentHighGrad ];
        startConnecetdIdx = [ startConnecetdIdx; ...        % update segment start boundaries
                currentHighGrad+1 ];
        idxPosGrad(1) = [];                                 % delete first selected index from the inspected list
        idxNegGrad(idxNegGrad<currentHighGrad) = [];        % delete smaller indexes from the other list
    elseif ~checkPosList  && ~isempty(idxNegGrad)  
        currentHighGrad = idxNegGrad(1);
        endConnecetdIdx = [ endConnecetdIdx ;...            % update segment end boundaries
              currentHighGrad ];
        startConnecetdIdx = [ startConnecetdIdx; ...     	% update segment start boundaries
                currentHighGrad+1 ];
        idxNegGrad(1) = [];                                 % delete first selected index from the inspected list
        idxPosGrad(idxPosGrad<currentHighGrad) = [];        % delete smaller indexes from the other list
    end
    checkPosList = ~checkPosList;
end

endConnecetdIdx = [endConnecetdIdx; length(xS)];
