function bin_mask = magicwandmono(im, ylist, xlist, tolerance)
% MAGICWAND simulates the Photoshop's magic wand tool
% It allows selection of connected pixels whose colors are
% within a defined tolerance of reference pixels.
%
% SYNTAX
%    bin_mask = magicwand(im, ylist, xlist, tolerance);
%
% INPUT
%   im:    input image RGB
%   ylist: vector of row cordinates    (reference pixels)
%   xlist: vector of column cordinates (reference pixels)
%   tolerance: distance to reference pixels
%
% OUTPUT
%   bin_mask: binary mask of selected regions
%
% EXAMPLE
%   The following code selects the girl's face:
%   im = imread('test.jpg');
%   bin_mask = magicwand(im, [199 217], [318 371], 50);
%   subplot(1, 2, 1); imshow(im);
%   subplot(1, 2, 2); imshow(bin_mask);
%
% NOTES
%   * Tested with MATLAB R13 & Image Processing Toolbox
%   * C++ version by Daniel Lau and Yoram Tal are also available in MATLAB Central
%
% (C) 2004 Son Lam Phung 
% Email: s.phung@ecu.edu.au 
% Edith Cowan University

H = size(im, 1); % image height
W = size(im, 2); % image width

% Check arguments
if any(ylist > H)
    s = sprintf('Row cordinates greater than the image height (%g).', H);
    error(s);
end

if any(xlist > W)
    s = sprintf('Column cordinates greater than the image height (%g).', W);
    error(s);
end

if ndims(im) ~=2
    error('Input image must be in BW format');
end

c_r = double(im(:, :, 1)); % Red channel
N = length(ylist); % Number of reference pixels

% Find all pixels whose colors fall within the specified tolerance
mono_mask = false(H, W);
for idx = 1:length(ylist)
    ref_r = double(im(ylist(idx), xlist(idx), 1));
    mono_mask = mono_mask | ...
                 ((c_r - ref_r) .^ 2) <= tolerance ^ 2;
end

% Connected component labelling
[objects, count] = bwlabel(mono_mask, 8); 

% Initialize output mask
bin_mask = false(H, W);

% Linear indices of reference pixels
pos_idxs = (xlist - 1) * H + ylist;

for idx = 1:count
    object = (objects == idx); % an object
   
    % Add to output mask if the object contains a reference pixel 
    if any(object(pos_idxs))
        bin_mask = bin_mask | object;
    end
end