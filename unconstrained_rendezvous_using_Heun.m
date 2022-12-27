%MATLAB implementation of 'Autonomous rendezvous using artificial potential
%function guidance' - unconstrained rendezvous usign Heun
%Edoardo Sampaolesi

clear; clc; close all;

%utils variables

tf = 15; % time
N = 300; %nodes
h = tf/N; %steps size

k = 0.35; %positive gain
p1 = 1;
p2 = 1;
p3 = 3;
P = diag([p1,p2,p3]);

%initial conditions
IV = [
    100 
    100 
    0 
    0.1 
    -0.002
    0];

%APF
V = @(x) x'*P*x;
%needed for drawing APF
V2 = @(x,y) p1*x.^2 + p2*y.^2;
%APF derivate
Vprimo = @(x,v) 2*x'*P*v;
%CW equation
CW = @(x,v) [
            3*x(1) + 2*v(2)
            -2*v(1)
            -x(3)
            ];

Pos = IV(1:3);
Vel = IV(4:6);
solV = V(Pos);
solVprimo = Vprimo(Pos,Vel);
time = [0 0];
index = 1;

for i = 1:N
     if solVprimo(i) < 0
           Vel(:,i+1) = Vel(:,i) + (h/2)*CW(Pos(:,i),Vel(:,i)) + (h/2)*CW(Pos(:,i),Vel(:,i) + h*CW(Pos(:,i),Vel(:,i)));
     else
          Vel(:,i+1) = -2*k*Pos(:,i)'*P;
     end
     Pos(:,i+1) = Pos(:,i) + (h/2)*Vel(:,i) + (h/2)*Vel(:,i+1);
     solV(i+1) = V(Pos(:, i+1)); %calc APF
     solVprimo(i+1) = Vprimo(Pos(:, i+1),Vel(:,i+1)); %calc APF derivated
     %calc when the chaser is near the taget
     if sqrt([1 1 1]*Pos(:,i).^2) < 1
         time(index) = i*h;
         index = 2;
     end
end


figure; set(gcf,'position',[10,10,1000,700])
set(0,'defaultTextInterpreter','latex');
%top sx plot
subplot(2,2,1);
plot3(Pos(1,:),Pos(2,:),Pos(3,:),'ko','MarkerSize',3); hold on; grid on;
plot3(Pos(1,1),Pos(2,1),Pos(3,1),'go','LineWidth',5);
plot3(0,0,0,'ro','LineWidth',3);
if Pos(3,1) == 0
    view(2);
end
legend('',sprintf('start (%i,%i,%i)',Pos(1,1),Pos(2,1),Pos(3,1)),'target','Location','best')
title('\textbf{Evolution 2d of $\xi$ , $\eta$ , $\zeta$ using Heun}',sprintf('Time: %i Nodes: %i Steps size: %0.5g',tf,N,h))
xlabel('$\xi$'); ylabel('$\eta$');zlabel('$\zeta$');
%top dx plot
subplot(2,2,2);
plot(0:h:tf,solV,'k-'); hold on;
plot(0:h:tf,solVprimo,'r--');
legend('V',"V'",'Location','best')
xlabel('Time'); ylabel("V and V'");
title("\textbf{Artificial potential function (black) and Rate of descent of potential (red)}",sprintf('Time to bring chaser to within 1 m to the target: %0.5g',time(1)))
%bottom sx plot
subplot(2,2,3);
plot(0:h:tf,Pos(1,:),'r-'); hold on;
plot(0:h:tf,Pos(2,:),'b-');
plot(0:h:tf,Pos(3,:),'k-');
legend('\xi','\eta','\zeta','location','best')
title('\textbf{Evolution of $\xi$ , $\eta$ , $\zeta$}')
xlabel('Time'); ylabel('Variables');
%bottom dx plot
subplot(2,2,4);
x = -100:0.1:100;  % define range and mesh of x and y which will be shown in figure
k = -100:0.1:100;
[X,Y] = meshgrid(x, k);
surf(X, Y, V2(X,Y),'EdgeColor','none'); hold on;
plot3(0,0,0,'ro','LineWidth',3); hold on; %origin
plot3(Pos(1,:),Pos(2,:),solV,'LineWidth',3,'Color','black');
view(-25,60);
legend('APF','target pos','V path','Location','best')
xlabel('$\xi$'); ylabel('$\eta$');zlabel('V');
title('\textbf{Unconstrained rendezvous}');