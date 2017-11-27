function OCTmask = getOCTmask(im)

%GETOCTMASK - returns a logical matrix of the same size of the first 2
%dimension of the input image and identifying the pixels not belonging to
%the scan.
%
%
% Author: Giovanni Ometto
% Work address: C274 Tait Building City, University of London, London, EC1V 0HB (UK)
% email: giovanni.ometto@city.ac.uk
% Website: http://www.city.ac.uk
% Jan 2017; Last revision:

% convert image to grayscale if not already
if size(im,3)>1, greyIm = rgb2gray(im); else, greyIm = im; end

coeff = size(greyIm,1)/372;

OCTmask = zeros(size(greyIm));

tolerance = 0.5;
se1 = strel('disk', round(2*coeff));
se2 = strel('disk', round(4*coeff));

% check which side has misalignment by looking at the variance (UPLEFT / UPRIGHT / LEFT / RIGHT)
LEFT = var(double(greyIm(:,1))) < 2;
RIGHT = var(double(greyIm(:,size(greyIm,2)))) < 2;
UPLEFT = var(double( greyIm(2, 1 : round(size(greyIm,2)*1/5))) ) < 2;
UPRIGHT = var(double( greyIm(2, round(size(greyIm,2)*4/5) : size(greyIm,2))) ) < 2;

coordsX = [] ;
coordsY = [] ; 

if LEFT||UPLEFT
    coordsX = [coordsX, 1];
    coordsY = [coordsY, 1];
end

if RIGHT||UPRIGHT
    coordsX = [coordsX, 1];
    coordsY = [coordsY, size(greyIm,2)];
end

% get the connected area
if ~isempty(coordsX)
    OCTmask = magicwandmono(greyIm, coordsX, coordsY, tolerance);
    OCTmask = imerode(OCTmask,se1);
    OCTmask = imdilate(OCTmask,se2);
    OCTmask = imerode(OCTmask,se1);
end

OCTmask = imfill(~OCTmask,'holes');

treshPseudoEnt = 0.008;

if sum(sum((gradient(OCTmask))))/(size(OCTmask,1)*size(OCTmask,2)) > treshPseudoEnt, OCTmask = repmat(var(double(greyIm),0,1)>2,size(greyIm,1),1); end


