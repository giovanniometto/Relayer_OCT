function [ILM] = snakeILM(im,ILMguess)

% Author: Giovanni Ometto
% Work address: C274 Tait Building City, University of London, London, EC1V 0HB (UK)
% email: giovanni.ometto@city.ac.uk
% Website: http://www.city.ac.uk
% Feb 2017; Last revision: Jun 2017

% convert image to grayscale if not already
if size(im,3)>1, greyIm = rgb2gray(im); else, greyIm = im; end

varianceTest = var(ILMguess(:,1));

% standardisation coefficient
% coeff = size(im,1)/372;
coeff =1;

% Set parameters and calculate the snake

Options.Wline=0; 
Options.Wedge=30; 
Options.Sigma1=round(3*coeff); 
Options.Sigma2=round(3*coeff);
Options.Alpha=.3; 
Options.Beta=100; 
Options.Gamma=1; 
Options.Kappa=.8;
Options.Iterations=round(50+(varianceTest)/100);
% Options.Verbose = true;
if Options.Iterations > 300, Options.Iterations = 300; end

ILM = Snake1D(greyIm,ILMguess,Options);
