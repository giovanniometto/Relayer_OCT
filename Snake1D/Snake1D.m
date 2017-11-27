function P = Snake1D(I,P,Options)
% Snake function limited to the vertical movement of the points.
%
% [O]=Snake1D(I,P,Options)
%  
% inputs,
%   I : An Image of type double
%   P : List with coordinates descriping the rough contour N x 2
%   Options : A struct with all snake options
%   
% outputs,
%   O : List with coordinates of the final contour M x 2
%
% options (general),
%  Option.Verbose : If true show important images, default false
%  Options.nPoints : Number of contour points, default pixel width of the image
%  Options.Gamma : Time step, default 1
%  Options.Iterations : Number of iterations, default 100
%
% options (Image Edge Energy / Image force))
%  Options.Sigma1 : Sigma used to calculate image derivatives, default 10
%  Options.Wline : Attraction to lines, if negative to black lines otherwise white
%                    lines , default 0.04
%  Options.Wedge : Attraction to edges, default 2.0
%  Options.Sigma2 : Sigma used to calculate the gradient of the edge energy
%                    image (which gives the image force), default 20
%
% options (Snake)
%  Options.Alpha : Membrame energy  (first order), default 0.2
%  Options.Beta : Thin plate energy (second order), default 0.2
%  Options.Kappa : Weight of external image force, default 2
%  
%   
% Options.Wline=0; Options.Wedge=20; Options.Sigma1=5; Options.Sigma2=8; Options.Alph=0; Options.Beta=0; Options.Gamma=1; Options.GIterations=0; Options.Iterations=100;
% [O,J]=Snake1D(I,P,Options);

% Process inputs
defaultoptions=struct('Verbose',false,'nPoints',size(P,1),'Wline',0.04,'Wedge',2,'Wterm',0.01,'Sigma1',10,'Sigma2',20,'Alpha',0,'Beta',0,'Delta',0.1,'Gamma',1,'Kappa',2,'Iterations',100,'GIterations',0,'Mu',0.2,'Sigma3',1);
if(~exist('Options','var'))
    Options=defaultoptions; 
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
         if(~isfield(Options,tags{i})), Options.(tags{i})=defaultoptions.(tags{i}); end
    end
    if(length(tags)~=length(fieldnames(Options)))
        warning('snake:unknownoption','unknown options found');
    end
end

% Convert input to double
I = im2double(I);

% If color image convert to grayscale
if(size(I,3)==3), I=rgb2gray(I); end

% Transform the Image into an External Energy Image
Eext = ExternalForceImage1D(I,Options.Wline, Options.Wedge, Options.Sigma1);

% Make the external force (flow) field.
Fy=ImageDerivatives1D(Eext,Options.Sigma2);
Fext(:,:,1)=-Fy*2*Options.Sigma2^2;

% Show the image, contour and force field
if(Options.Verbose)
    h=figure(1); set(h,'render','opengl')
     subplot(2,2,1),
      imshow(I,[]); 
      hold on; plot(P(:,2),P(:,1),'b.'); hold on;
      title('The image with initial contour')
     subplot(2,2,2),
      imshow(Eext,[]); 
      title('The external energy');
     subplot(2,2,3), 
      [x,y]=ndgrid(1:5:size(Fext,1),1:5:size(Fext,2));
      FextX = Fext(1:5:end,1:5:end,1);
      imshow(I), hold on; quiver(y,x,zeros(size(FextX)),FextX);
      title('The external force field ')
     subplot(2,2,4), 
      imshow(I), hold on; plot(P(:,2),P(:,1),'b.'); 
      title('Snake movement ')
    drawnow
end


% Make the interal force matrix, which constrains the moving points to a
% smooth contour
S=SnakeInternalForceMatrix1D(Options.nPoints,Options.Alpha,Options.Beta,Options.Gamma);
h=[];
for i=1:Options.Iterations
    P=SnakeMoveIteration1D(S,P,Fext,Options.Gamma,Options.Kappa);

    % Show current contour
    if(Options.Verbose)
        if(ishandle(h)), delete(h), end
        h=plot(P(:,2),P(:,1),'r.');
        c=i/Options.Iterations;
        plot(P(:,2),P(:,1),'-','Color',[c 1-c 0]);  drawnow
    end
end

P = flipud(P');
