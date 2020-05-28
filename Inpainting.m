clear;

% Choose portion of remaining pixels:
q = 0.03;

% Import image as grayscale:
I = rgb2gray(imread('Lenna.png'));
% imshow(I)
% pause

% Choose parameters:
nx = size(I,1);
ny = size(I,2);
StopTime = 500;
dx = 1;
dy = 1;
dt = (dx^2+dy^2)/8;
L = dt/dx^2;

% Randomly delete (1-q) of the image;
r = unique(sort(randi([1 nx*ny],1,round(q*nx*ny))));

N = (zeros(size(I)));
N(r)=I(r);
Percent_of_Pixels = 100*nnz(N)/numel(N);
disp(Percent_of_Pixels);
figure
imshow(uint8(N))
pause

%%

% Run diffusion:
U = N;
V = ones(size(N));
Z = zeros(nx,ny);

while max(max(abs(V-U)))/max(max(abs(U)))>1e-3
    V=U;
    U1 = circshift(U,1,1);
    U2 = circshift(U,1,2);
    U1inv = circshift(U,-1,1);
    U2inv = circshift(U,-1,2);
    
    % Equation for U:
    U = (1-4*L)*U + L*(U2inv + U2 + U1inv + U1) ;

    %Neumann Boundary Conditions:
    U(:,1) = U(:,3);     
    U(:,ny) = U(:,ny-2);   
    U(1,:) = U(3,:);     
    U(nx,:) = U(nx-2,:);   
    
    % "Boundary condition" on the known pixels:
    U(N>1) = N(N>1);

    frame = (uint8(U));
    imshow(uint8(U))
end
