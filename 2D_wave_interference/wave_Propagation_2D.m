%% 2d
close all
clear all
clc,
L=1000;% L=2000 metres
T=0.5; % T=4 sec 
c=2500; % c=1500 m//s
f0=25; % f0=5 HZ 

lamda_Min=c/(2.5*f0);
dx= 10;
dy=10;
dt=dx/c/sqrt(2); % Pas d echantillonage
kx=round(L/dx)+1;
ky=round(L/dy)+1; % Meme pas d'echantilonage pour x et y 
kt=round(T/dt)+1;

t0=1.5*(6^0.5)/(pi*f0);
t=0:dt:T;
r=100*(1-2.*(pi.^2).*(f0.^2).*(t-t0).^2).*exp((-pi.^2).*(f0.^2).*((t-t0).^2)); 

p=zeros(kt,kx,ky);
Mx=round(L/dx)+1;
My=Mx;
x=linspace(0,L,Mx); % Distance along x
y=linspace(0,L,My); % Distance along y
sigma= ones(1,kx);
 width=20;
 for kk=1:width
      sigma(kk) = exp(-0.3*((width-kk)/width)^2);
      sigma(kx-kk+1)=sigma(kk);
 end
 gamma=sigma'*sigma;
 
%  figure(1) 
% plot(sigma);
% 
% figure(2)
% imagesc(gamma);

 p(1,:,:)=0;
 p(2,:,:)=0;
 
xs=500;
ys=250;
is=round(xs/dx)+1;
js= round(ys/dy)+1;

x1s=500;
y1s=750;
i1s=round(x1s/dx)+1;
j1s= round(y1s/dy)+1;

 for i=2:kt-1
      for j=2:kx-1
        for d=2:ky-1
        p(i+1,j,d)=(((c^2*dt^2)/(dx^2))*(p(i,j+1,d)+p(i,j,d+1)-4*p(i,j,d)+p(i,j-1,d)+p(i,j,d-1))+2*p(i,j,d)-p(i-1,j,d));
        end
      end
      p(i+1,is,js)=p(i+1,is,js)+dt^2*r(i);
      p(i+1,i1s,j1s)=p(i+1,i1s,j1s)+dt^2*r(i);

  for j= 1:kx
    for d= 1:ky
    p(i+1,j,d)= p(i+1,j,d)*gamma(j,d);
    p(i,j,d)= p(i,j,d)*gamma(j,d);
    end
  end  
 end
    p(:,j,1)=0;
    p(:,j,end)=0;  
    p(:,1,d)=0;
    p(:,end,d)=0;   

myVideo = VideoWriter('myVideoFile'); %open video file
myVideo.FrameRate = 5;  %can adjust this, 5 - 10 works well for me
open(myVideo)

 figure(3)
      
 for i=1:kt-1
      subplot(211)
     mesh(x,y,squeeze(p(i,:,:))); colormap(flipud(bone)); colorbar; caxis([-0.0002 0.0002]); axis([0 L 0 L -0.001 0.001]);
     title(sprintf('wave interference simulation _ 3D view _ (Time=%.2f)',t(1,i))); 
     pause(0.01);
 
      subplot(212)
     sc=reshape(p(i,:,:),kx,ky);
     imagesc(x,y,sc); colormap(flipud(bone)); colorbar; caxis([-0.0002,0.0002]);
     title(sprintf('wave interference simulation 2D view _ (Time=%.2f)',t(1,i))); 
     pause(0.01);
 frame = getframe(gcf); %get frame
   writeVideo(myVideo, frame);
end
close(myVideo) 
 
 
 
 




