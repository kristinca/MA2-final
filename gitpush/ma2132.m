%Naravna konvekcija

clc
clear 

nx=40; ny=40; h=1/nx;
dt=0.00005;nstep=200000;
tlakN=180;PP=1.8;
u=zeros(nx+1,ny+2); v=zeros(nx+2,ny+1); p=zeros(nx+2,ny+2);
ut=zeros(nx+1,ny+2); vt=zeros(nx+2,ny+1); c=zeros(nx+2,ny+2)+0.25;
c(2,3:ny)=1/3;c(nx+1,3:ny)=1/3;c(3:nx,2)=1/3;c(3:nx,ny+1)=1/3;
c(2,2)=1/2;c(2,ny+1)=1/2;c(nx+1,2)=1/2;c(nx+1,ny+1)=1/2;

T0=20;
T=T0*zeros(nx+2,ny+2); nT=T0*zeros(nx+2,ny+2);
% lahko tudi T0*ones(nx+2,ny+2) za doloceno temperaturo T0 povsod po domenu
% parametri variant sem podala v Excel tabeli z imenom m
% to je prikaz dostopa do teh parametrov za varianto A
% pri tem m importiramo v MATLAB kot numeric matrix
Alfa1=m(1,1); Alfa2=m(1,2);
Beta1=m(1,3); Beta2=m(1,4);
mu1=m(1,5); mu2=m(1,6);

for is=1:nstep
%rp, hitrost je enaka nic po robovih
    u(1:nx+1,1)= zeros(nx+1,1);
    u(1:nx+1,ny+2)=zeros(nx+1,1);
    v(1,1:ny+1)= zeros(1,ny+1);
    v(nx+2,1:ny+1)=zeros(1,ny+1);
    nT(1:nx+2,1)=20; %levi rob segrevamo
    nT(1:nx+2,ny+2)=zeros; %desni rob ima temperaturo nic 
        for i=2:nx
            for j=2:ny+1
                %rešujemo za hitrost ob tem trenutku
                if i<=nx/2
                    ut(i,j)=u(i,j)+dt*(-(0.25/h)*((u(i+1,j)+u(i,j))^2-(u(i,j)+...
                    u(i-1,j))^2+(u(i,j+1)+u(i,j))*(v(i+1,j)+...
                    v(i,j))-(u(i,j)+u(i,j-1))*(v(i+1,j-1)+v(i,j-1)))+...
                    (mu1/h^2)*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4*u(i,j))+Beta1*9.81*(T(i,j)-T0));
                else
                    ut(i,j)=u(i,j)+dt*(-(0.25/h)*((u(i+1,j)+u(i,j))^2-(u(i,j)+...
                    u(i-1,j))^2+(u(i,j+1)+u(i,j))*(v(i+1,j)+...
                    v(i,j))-(u(i,j)+u(i,j-1))*(v(i+1,j-1)+v(i,j-1)))+...
                    (mu2/h^2)*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4*u(i,j))+Beta2*9.81*(T(i,j)-T0));
                end
            end
        end
        for i=2:nx+1
            for j=2:ny
            %rešujemo za hitrost ob tem trenutku
                if i<=nx/2
                    vt(i,j)=v(i,j)+dt*(-(0.25/h)*((u(i,j+1)+u(i,j))*(v(i+1,j)+...
                    v(i,j))-(u(i-1,j+1)+u(i-1,j))*(v(i,j)+v(i-1,j))+...
                    (v(i,j+1)+v(i,j))^2-(v(i,j)+v(i,j-1))^2)+...
                    (mu1/h^2)*(v(i+1,j)+v(i-1,j)+v(i,j+1)+v(i,j-1)-4*v(i,j)));
                else
                    vt(i,j)=v(i,j)+dt*(-(0.25/h)*((u(i,j+1)+u(i,j))*(v(i+1,j)+...
                    v(i,j))-(u(i-1,j+1)+u(i-1,j))*(v(i,j)+v(i-1,j))+...
                    (v(i,j+1)+v(i,j))^2-(v(i,j)+v(i,j-1))^2)+...
                    (mu2/h^2)*(v(i+1,j)+v(i-1,j)+v(i,j+1)+v(i,j-1)-4*v(i,j)));
                end
            end
        end
        for ii=1:tlakN
        %rešujemo za tlaka
            for i=2:nx+1
                for j=2:ny+1
                    p(i,j)=PP*c(i,j)*(p(i+1,j)+p(i-1,j)+p(i,j+1)+p(i,j-1)-...
                    (h/dt)*(ut(i,j)-ut(i-1,j)+vt(i,j)-vt(i,j-1)))+(1-PP)*p(i,j);
                end
            end
        end
        %rešitev za tlaka uporabimo za dolocitev hitrosti
        u(2:nx,2:ny+1)=...
        ut(2:nx,2:ny+1)-(dt/h)*(p(3:nx+1,2:ny+1)-p(2:nx,2:ny+1));
        v(2:nx+1,2:ny)=...
        vt(2:nx+1,2:ny)-(dt/h)*(p(2:nx+1,3:ny+1)-p(2:nx+1,2:ny));
           for i=2:nx+1
                for j=2:ny+1
                    if i<=(nx+2)/2
                        nT(i,j)=...
                            T(i,j)+...
                            dt*(...
                            (Alfa2/h^2)*(T(i,j+1)+T(i,j-1)+T(i+1,j)+T(i-1,j)-4*T(i,j))-...
                            (1/(2*h))*(...
                             v(i,j)*(T(i,j+1)-T(i,j-1))+...
                             u(i,j)*(T(i+1,j)-T(i-1,j))...
                             ));
                    end
                    nT(1,2:ny+2)=nT(2,2:ny+2); nT(nx+2,2:ny+2)=nT(nx+1,2:ny+2);
                    t=t+dt;
                    T=nT;
                end
           end
end