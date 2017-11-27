function Fext=GVFOptimizeImageForces1D(Fext, Mu, Iterations, Sigma)
% Gradient vector flow (GVF) on the vertical direction
%
% Fext = GVFOptimizeImageForces1D(Fext, Mu, Iterations, Sigma) 
% 
% inputs,
%   Fext : The image force vector field N x M x 2
%   Mu : Is a trade of scalar between noise and real edge forces
%   Iterations : The number of GVF itterations
%   Sigma : Used when calculating the Laplacian
% 
% outputs,
%   Fext : The GVF optimized image force vector field
%
% Author: Giovanni Ometto
% Work address: C274 Tait Building City, University of London, London, EC1V 0HB (UK)
% email: giovanni.ometto@city.ac.uk
% Website: http://www.city.ac.uk
% Feb 2017; Last revision:

% Squared magnitude of force field
Fy= Fext(:,:,1);

% Calculate magnitude
sMag = Fy.^2;

% Set new vector-field to initial field
u=Fy;
  
% Iteratively perform the Gradient Vector Flow (GVF)
for i=1:Iterations,
  % Calculate Laplacian
  Uyy=ImageDerivatives2D(u,Sigma,'yy');
  u = u + Mu*(Uyy) - sMag.*(u-Fy);
  
end

Fext(:,:,1) = u;
