function J=ImageDerivatives1D(I,sigma)
% Gaussian based image derivatives
%
%  J=ImageDerivatives1D(I,sigma)
%
% inputs,
%   I : The image
%   sigma : Gaussian Sigma
%
% outputs,
%   J : The image vertical derivative
%
% Author: Giovanni Ometto
% Work address: C274 Tait Building City, University of London, London, EC1V 0HB (UK)
% email: giovanni.ometto@city.ac.uk
% Website: http://www.city.ac.uk
% Feb 2017; Last revision:


% Make derivatives kernels
[y,x]=ndgrid(floor(-3*sigma):ceil(3*sigma),floor(-3*sigma):ceil(3*sigma));
DGauss=-(y./(2*pi*sigma^4)).*exp(-(x.^2+y.^2)/(2*sigma^2));
J = imfilter(I,DGauss,'conv','symmetric');