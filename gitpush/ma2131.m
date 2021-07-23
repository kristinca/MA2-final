%Prenos toplote s prevodom

clc
clear 

nx=30;ny=30; dt=0.0005; nstep=2000; h=1/nx;
T=zeros(nx+2,ny+2); k=0.00019; 
time=0;

for is=1:nstep
%rp
    nT(1:nx+2,1)=20; %lev rob T=20
    nT(1:nx+2,ny+2)=zeros; %desni rob T=0
        for i=2:nx+1
            for j=2:ny+1
                nT(i,j)=T(i,j)+k*(1/h^2)*(T(i,j+1)+T(i,j-1)+T(i+1,j)+T(i-1,j)-4*T(i,j));
            end
        end
        nT(1,2:ny+1)=nT(2,2:ny+1); %ni prenosa toplote skozi roba zgori
        nT(nx+2,2:ny+1)=nT(nx,2:ny+1); %ni prenosa toplote skozi roba spodi
        T=nT; % kot vhodni podatek za temperaturo dolocimo trenutno vrednost temperature
        time=time+dt;
end

% plot
figure (1)
P=imagesc(T);
title(sprintf('Prenos toplote s prevodom, dt = %0.4f',dt));
xlabel('x'); ylabel('y');
colorbar
% set(gca,’ytick’,[]); set(gca,’xtick’,[]); %brez oznake o elementih na x in y osi