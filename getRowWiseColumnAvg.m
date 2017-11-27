function verticalProfile = getRowWiseColumnAvg(bwim,xCoord,halfwidth)

%GETVERTICALTILEPROFILE - Returns an array corresponding to the averaged
%values of a BW image computed along the rows of a tile described by the
%given centre, width (+1), height (+1) and angle (in radians)
%
%
% Author: Giovanni Ometto
% Work address: C274 Tait Building City, University of London, London, EC1V 0HB (UK)
% email: giovanni.ometto@city.ac.uk
% Website: http://www.city.ac.uk
% Jan 2017; Last revision:

% check that it is a bw image and if not make it one
if size(bwim,3) ~= 1
    bwim = rgb2gray ( bwim );
end

halfTileW = halfwidth; 
imw = size(bwim,2);
verticalProfile = mean(bwim(:,max(1,xCoord-halfTileW): min(imw,xCoord+halfTileW)), 2);
  
end

