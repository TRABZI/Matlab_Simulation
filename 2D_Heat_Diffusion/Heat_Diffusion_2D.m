close all
clear all
clc,

% 2D Heat Diffusion in metalic plate
rho=1; % density 
%k= 2; % Thermal conductivity 
cp=5; % specific heat 
%A=k/rho/cp;

% Descritisation of the metal plate
Lx=10;
Ly=10;
Nx=101; Nt=500;
Ny=101;
dx=Lx/(Nx-1); % step length
dy=Ly/(Ny-1);

% Time step 
c=1;
C=0.05;
dt=C*dx/c; % dt 

% Field variables 
Tn=zeros(Ny,Nx);  % Temperature
x=linspace(0,Lx,Nx); % Distance along x
y=linspace(0,Ly,Ny); % Distance along y
[X,Y]=meshgrid(x,y);
k=ones(Ny,Nx);
% k(20:25, 30:35)=0.000001; 
% figure(1)
% imagesc(k)
% Defining the initiale conditions 
Tn(:,:)=0;
t=0 ;
% loop 

myVideo = VideoWriter('VideoFile'); %open video file
myVideo.FrameRate = 5;  %can adjust this, 5 - 10 works well for me
open(myVideo)

figure(2)

for n=1:Nt
    T=Tn; % Save the temperature into T for using it later 
    t=t+dt; 
    % loop for the new temperature for defferent x positions
    for i=2:Nx-1
        for j=2:Ny-1
         Tn(j,i)=T(j,i)+...
             dt*(k(j,i)/rho/cp)*...
             ((T(j,i+1)+T(j+1,i)-4*T(j,i)+T(j,i-1)+T(j-1,i))/dx/dx);
        end
    end 
    %BC:
    Tn(1,:)=0; Tn(end,:)=0; % Dirichlet condition 
    Tn(:,1)=100; Tn(:,end)=T(:,end-1);
     
%plot
% subplot(211)
imagesc(Tn); colormap;
title(sprintf('2D Heat diffusion in metal plate (Time=%f seconds)',t)); 
colorbar;
% subplot(212)
% mesh(x,y,Tn), axis([0 Lx 0 Ly 0 150]);
% colorbar;
% title(sprintf('Time=%f seconds',t)); 
pause(0.05);
frame = getframe(gcf); %get frame
   writeVideo(myVideo, frame);
end 
close(myVideo) 

