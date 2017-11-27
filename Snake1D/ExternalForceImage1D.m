function Eextern = ExternalForceImage1D(I, Wline, Wedge, Sigma)
% Calculate the y external forces (termination functional does not apply
% to this case)
% 
% inputs, 
%  I : The image
%  Sigma : Sigma used to calculated image derivatives 
%  Wline : Attraction to lines, if negative to black lines otherwise white
%          lines
%  Wedge : Attraction to edges
%
% outputs,
%  Eextern : The energy function described by the image
%
% Author: Giovanni Ometto
% Work address: C274 Tait Building City, University of London, London, EC1V 0HB (UK)
% email: giovanni.ometto@city.ac.uk
% Website: http://www.city.ac.uk
% Feb 2017; Last revision:

Iy=ImageDerivatives1D(I,Sigma);

Eline = imgaussian(I,Sigma);
Eedge = sqrt(Iy.^2); 

Eextern= (Wline*Eline - Wedge*Eedge ); 

