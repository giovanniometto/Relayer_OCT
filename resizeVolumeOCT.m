function resizedVolume = resizeVolumeOCT(octVolume, factor, dimension)

% Author: Giovanni Ometto
% Work address: C274 Tait Building City, University of London, London, EC1V 0HB (UK)
% email: giovanni.ometto@city.ac.uk
% Website: http://www.city.ac.uk
% Jan 2017; Last revision: Oct 2017

if nargin < 3
    dimension = 1;
end

if dimension == 1
    [X1, Y1, Z1]  = meshgrid (1:size(octVolume,2), 1:size(octVolume,1), 1:size(octVolume,3));
    [X2, Y2, Z2]  = meshgrid (1:size(octVolume,2), 1:1/factor:size(octVolume,1), 1:size(octVolume,3));
    resizedVolume = uint8(round(interp3( X1, Y1, Z1, double(octVolume), X2, Y2, Z2)));
elseif dimension == 2
    [X1, Y1, Z1]  = meshgrid (1:size(octVolume,2), 1:size(octVolume,1), 1:size(octVolume,3));
    [X2, Y2, Z2]  = meshgrid (1:1/factor:size(octVolume,2), 1:size(octVolume,1), 1:size(octVolume,3));
    resizedVolume = uint8(round(interp3( X1, Y1, Z1, double(octVolume), X2, Y2, Z2)));
else 
    return
end


