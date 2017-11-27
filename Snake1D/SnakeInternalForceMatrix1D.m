function B=SnakeInternalForceMatrix1D(nPoints,alpha,beta,gamma)
%
% B=SnakeInternalForceMatrix2D(nPoints,alpha,beta,gamma)
%
% inputs,
%   nPoints : The number of snake contour points
%   alpha : membrame energy  (first order)
%   beta : thin plate energy (second order)
%   gamma : Step Size (Time)
%
% outputs,
%   B : The Snake Smoothness regulation matrix
%
% Function is written by D.Kroon University of Twente (July 2010)

% Penta diagonal matrix, one row:
b(1)=beta;
b(2)=-(alpha + 4*beta);
b(3)=(2*alpha + 6 *beta);
b(4)=b(2);
b(5)=b(1);

diagNPoints = eye(nPoints);

% Make the penta matrix (for every contour point)
A=b(1)*circshift(diagNPoints,2);
A=A+b(2)*circshift(diagNPoints,1);
A=A+b(3)*circshift(diagNPoints,0);
A=A+b(4)*circshift(diagNPoints,-1);
A=A+b(5)*circshift(diagNPoints,-2);

% Calculate the inverse
B=inv(A + gamma.* eye(nPoints));

% Modify matrix adding the values in top right and bottom left corners to
% the top left and bottom right corners respectively of matrix B - to get
% the snake of an open shape
T = ones(size(B));
B = B + fliplr(triu(B,round(nPoints/2)));   % add top right to top left
B(logical(triu(T,round(nPoints/2)))) = 0;   % remove top right values
B = B + fliplr(tril(B,-round(nPoints/2)));  % add bottom left to bottom right
B(logical(tril(T,-round(nPoints/2)))) = 0;  % remove bottom left values
