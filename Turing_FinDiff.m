clear;

J = 80;
K = 11*J^2;
T = 1/10;
StopTime = 10*K;
dt = T/K;
dx = 1/J;
L = dt/dx^2;

a = -0.55;
b = 1.9;
gamma = 600;
d = 4.8;

[XX,YY] = meshgrid(dx*(1:J), dx*(1:J));
U = 1.35*ones(J,J) + rand(J,J)*1e-2;
V = 1.0425*ones(J,J)+ rand(J,J)*1e-2;

% v = VideoWriter('TuringFinDiffx.avi');
% open(v);
for k=1:StopTime-1
    U1 = circshift(U,1,1);
    U2 = circshift(U,1,2);
    U1inv = circshift(U,-1,1);
    U2inv = circshift(U,-1,2);
    
    % Equation for U:
    U = (1-4*L)*U + L*U2inv + L*U2 + L*U1inv + L*U1... 
                + gamma*dt*(a - U + U.^2.*V) ;

    U(:,1) = U(:,3);     %Neumann Boundary Conditions
    U(:,J) = U(:,J-2);   
    U(1,:) = U(3,:);     
    U(J,:) = U(J-2,:);   
    
    V1 = circshift(V,1,1);
    V2 = circshift(V,1,2);
    V1inv = circshift(V,-1,1);
    V2inv = circshift(V,-1,2);
    
    % Equation for V:
    V = (1-4*L*d)*V + L*d*V2inv + L*d*V2 + L*d*V1inv + L*d*V1... 
                + gamma*dt*(b - U.^2.*V) ;

    V(:,1) = V(:,3);     %Neumann Boundary Conditions
    V(:,J) = V(:,J-2);   
    V(1,:) = V(3,:);     
    V(J,:) = V(J-2,:);   
    
    if k/J^2==ceil(k/J^2)
        imagesc(U)
        caxis([0 2.2]);
        colorbar;
        shading interp;
        axis off
        drawnow;
        StopTime-1-k
    %     frame = getframe(gcf);
    %     writeVideo(v,frame);
    end
end
% close(v);
