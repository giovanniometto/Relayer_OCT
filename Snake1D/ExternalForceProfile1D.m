function Eextern = ExternalForceProfile1D(P,Wline, Wedge, Sigma)
% Eextern = ExternalForceProfile1D(I,Wline, Wedge, Wterm,Sigma)
% 
% inputs, 
%  P : The profile
%  Sigma : Sigma used to calculated image derivatives 
%  Wline : Attraction to peaks, if negative to negative otherwise positive
%          peaks
%  Wedge : Attraction to gradients
%
% outputs,
%  Eextern : The energy function described by the profile
%
% Author: Giovanni Ometto
% Work address: C274 Tait Building City, University of London, London, EC1V 0HB (UK)
% email: giovanni.ometto@city.ac.uk
% Website: http://www.city.ac.uk
% Jan 2017; Last revision:

Eline = imgaussian(I,Sigma);
Eedge = gradient(P); 

Eextern= (Wline*Eline - Wedge*Eedge); 

