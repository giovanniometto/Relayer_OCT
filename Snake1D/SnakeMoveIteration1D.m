function P=SnakeMoveIteration1D(B,P,Fext,gamma,kappa)
% One iteration of contour Snake vertical movement
%
% P=SnakeMoveIteration1D(S,P,Fext,gamma,kappa)
%
% inputs,
%   B : Internal force (smoothness) matrix
%   P : The contour points N x 2;
%   Fext : External vector field (from image)
%   gamma : Time step
%   kappa : External (image) field weight
%
% outputs,
%   P : The (vertically moved) contour points
%
% Author: Giovanni Ometto
% Work address: C274 Tait Building City, University of London, London, EC1V 0HB (UK)
% email: giovanni.ometto@city.ac.uk
% Website: http://www.city.ac.uk
% Feb 2017; Last revision:

% Clamp contour to boundary
P(:,1)=min(max(P(:,1),1),size(Fext,1));

% Get image force on the contour points
Fext1(:,1)=kappa*interp2(Fext(:,:,1),P(:,2),P(:,1));

% Nans to zeros
Fext1(isnan(Fext1))=0;

% generate Weighting vector for Fext1
% maxX = 20;
% x = 0 : (maxX) / (size(P,1)-1) : (maxX);
% heightWV = 5;
% weightingV = atan(x)';
% halfWeightingV =  weightingV(1:round(size(P,1)/2));
% weightingV(end-size(halfWeightingV)+1:end) = flipud(halfWeightingV);
% weightingV = 1 + ( (heightWV-1) * ( 1 - weightingV./max(weightingV) ));

% sizeCorrection = round(size(P,1)/30);
% x = pi/2:pi/2/(sizeCorrection-1):pi;
% heightCos = 10;
% weightingV = ones(size(P,1),1); 
% cosV = 1+heightCos+heightCos*cos(x);
% weightingV(1:sizeCorrection,1) = cosV;
% weightingV(end-sizeCorrection+1:end) = fliplr(cosV);

% Update
% ssx = gamma*P(:,1) + weightingV.*Fext1(:,1);
ssx = gamma*P(:,1) + Fext1(:,1);
P(:,1) = B * ssx;

% Clamp contour to boundary
P(:,1)=min(max(P(:,1),1),size(Fext,1));
    
