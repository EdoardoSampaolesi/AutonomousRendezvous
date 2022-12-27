%MATLAB implementation of 'Autonomous rendezvous using artificial potential
%function guidance' - terminal guidance rendezvous with obstacles
%Edoardo Sampaolesi

clear; clc; close all;

%utils variables

tf = 30; % time
N = 300; %nodes
h = tf/N; %steps size

k = 0.35; %positive gain
p1 = 1;
p2 = 1;
p3 = 3;
P = diag([p1,p2,p3]);

m1 = 0.05;
m2 = 0.05;
m3 = 0.05;
M = diag([m1,m2,m3]);

a1 = 30;
b1 = 21;
b = 40;
a = a1/2;

%initial conditions
IV = [
    100
    100 
    0 
    0.1 
    -0.002
    0];

PosObstacles = [ %1 %2 %3
                 50 2 80 %csi
                 40 50 75 %eta
                 0   0   0 %zeta
                ];
parameters = [ %1  %2   %3
               25000 25000 25000 %psi
                5      5     5  %sigma
              ];

size = length(PosObstacles(1,:));

o = 0:0.001:pi;
d1 = sqrt((a^2).*(cos(o) - 1).^2 + (b1^2).*sin(o).^2);
o = 2*pi:0.001:3*pi;
d2 = sqrt((a^2).*(cos(o) + 1).^2 + (b1^2).*sin(o).^2);
o = pi:0.001:2*pi;
d3 = sqrt( (a1^2).*cos(o).^2 + 4*(b^2).*sin(o).^2 );

r1 = max([d1 d2 d3]);

LAMBDA = @(i,Pos) parameters(1,i)*exp(-(parameters(2,i)^-1)*(Pos-PosObstacles(:,i))'*M*(Pos-PosObstacles(:,i)));

%APF
VUnconstrained = @(x) x'*P*x;
%APF derivate
VprimoUnconstrained = @(x,v) 2*x'*P*v;
VprimoObstacles = @(x,v) (2*x'*P - lambda(M,x,LAMBDA,size,parameters,PosObstacles))*v;
%CW equation
CW = @(x,v) [
            3*x(1) + 2*v(2)
            -2*v(1)
            -x(3)
            ];

%for drawing APF
xsurf = -100:0.1:100;  % define range and mesh of x and y which will be shown in figure
ysurf = -100:0.1:100;
[X,Y] = meshgrid(xsurf, ysurf);
V2 = p1*X.^2 + p2*Y.^2;

Pos = IV(1:3);
Vel = IV(4:6);
solV = VObstacles(Pos,LAMBDA,size,P);
solVprimo = VprimoObstacles(Pos,Vel);
time = [0 0];
index = 1;

clear o;
for i = 1:N
    distance = sqrt(Pos(1,i)^2 + Pos(2,i)^2 );
    if distance < r1
         alfa = atan( Pos(2,i) / Pos(1,i) ); 
         if(Pos(2,i) < 0)
             d = @(o)  sqrt( (a1^2)*cos(o)^2 + 4*(b^2)*sin(o)^2 );
             o = atan(Pos(2,i)/Pos(1,i));
             if(Pos(1,i) > 0)
                 tau = [
                       -a1*sin(o)
                        b*cos(o)
                        0
                        ];
             else
                 tau = [
                       a1*sin(o)
                        -b*cos(o)
                        0
                        ];
             end
         else 
             if(Pos(1,i) > 0)
                 o = atan( Pos(2,i) / abs(Pos(1,i)+a) ); 
                 tau = [
                         -a*sin(o)
                         b1*cos(o)
                         0
                    ];
                d = @(o)  sqrt((a^2)*(cos(o) + 1)^2 + (b1^2)*sin(o)^2);
             else
                 o = atan( Pos(2,i) / Pos(1,i) + (3/2)*pi); 
                 tau = [
                        a*sin(o)
                        -b1*cos(o)
                        0
                    ];
                d = @(o) sqrt((a^2)*(cos(-o) + 1)^2 + (b1^2)*sin(-o)^2);
             end
         end
         grad = Pos(:,i)'*P;
         if( distance < d(alfa)*1.1 )
            Vel(:,i+1) = k*tau * 2*sqrt(grad(1)^2 + grad(2)^2)/sqrt(tau(1)^2 + tau(2)^2);
         else
             if solVprimo(i) < 0
                   Vel(:,i+1) = Vel(:,i) + (h/2)*CW(Pos(:,i),Vel(:,i)) + (h/2)*CW(Pos(:,i),Vel(:,i) + h*CW(Pos(:,i),Vel(:,i)));
             else
                  Vel(:,i+1) = -2*k*grad;
             end
         end
    else
        if solVprimo(i) < 0
           Vel(:,i+1) = Vel(:,i) + (h/2)*CW(Pos(:,i),Vel(:,i)) + (h/2)*CW(Pos(:,i),Vel(:,i) + h*CW(Pos(:,i),Vel(:,i)));
         else
              Vel(:,i+1) = -k*(2*Pos(:,i)'*P - lambda(M,Pos(:,i),LAMBDA,size,parameters,PosObstacles));
         end
    end

     Pos(:,i+1) = Pos(:,i) + (h/2)*Vel(:,i) + (h/2)*Vel(:,i+1);

     if distance < r1
        solV(i+1) = VUnconstrained(Pos(:, i+1)); %calc APF
        solVprimo(i+1) = VprimoUnconstrained(Pos(:, i+1),Vel(:,i+1)); %calc APF derivated
     else
         solV(i+1) = VObstacles(Pos(:, i+1),LAMBDA,size,P); %calc APF
         solVprimo(i+1) = VprimoObstacles(Pos(:, i+1),Vel(:,i+1)); %calc APF derivated
     end
     
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

t = linspace(0,pi) ;
x = a*cos(t) - a;
y = b1*sin(t) ;
plot(x,y,'r'); hold on;
axis equal;

t = linspace(0,pi) ;
x = a*cos(t) + a;
y = b1*sin(t);
plot(x,y,'r'); hold on;

t = linspace(pi,2*pi) ;
x = a1*cos(t);
y = b*sin(t);
plot(x,y,'r'); hold on;
for i = 1:size
    plot3(PosObstacles(1,i),PosObstacles(2,i),PosObstacles(3,i),'bo','LineWidth',3); hold on;
end
if Pos(3,1) == 0
    view(2);
end
legend('',sprintf('start (%i,%i,%i)',Pos(1,1),Pos(2,1),Pos(3,1)),'target','constraint boundary','','','obstacles','Location','best')
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
surf(X, Y, V2,'EdgeColor','none'); hold on;
plot3(0,0,0,'ro','LineWidth',3); hold on; %origin
plot3(Pos(1,:),Pos(2,:),solV,'LineWidth',3,'Color','black'); hold on;
for i = 1:size
    gauss = parameters(1,i) * exp( -(parameters(2,i)^-1) * ( m1*(X-PosObstacles(1,i)).^2 + m2*(Y-PosObstacles(2,i)).^2 ) );
    surf(X, Y, gauss,'EdgeColor','none'); hold on;
end
view(-25,60);
legend('APF','target pos','V path','Location','best')
xlabel('$\xi$'); ylabel('$\eta$');zlabel('V');
title('\textbf{Terminal guidance rendezvous}');

%lambda values
function lmbd = lambda(M,Pos,LAMBDA,size,parameters,PosObstacles) 
    lmbd = [0 0 0];
    for i = 1:size
        lmbd(1) = lmbd(1) + M(1,1)*LAMBDA(i,Pos)*(parameters(2,i)^-1)*(Pos(1)-PosObstacles(1,i)); 
        lmbd(2) = lmbd(2) + M(2,2)*LAMBDA(i,Pos)*(parameters(2,i)^-1)*(Pos(2)-PosObstacles(2,i));
        lmbd(3) = lmbd(3) + M(3,3)*LAMBDA(i,Pos)*(parameters(2,i)^-1)*(Pos(3)-PosObstacles(3,i)); 
    end
end

%APF
function v = VObstacles(x,LAMBDA,size,P)
    v = x'*P*x;
    for i = 1:size
        v = v + LAMBDA(i,x);
    end
end