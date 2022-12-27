%MATLAB implementation of Clohessy-Wiltshire equations
%using fEul, Heun and solutions
%Edoardo Sampaolesi

clear; clc; close all;

%utils variables
n = 0.00113; % il low orbit n has this value

tf = 8000; % time
N = 80; %nodes
h = tf/N; %steps size

%soluzioni in forma chiusa
%solved CV equations
M = @(t) [
                  4-3*cos(n*t)      0      0            
                6*(sin(n*t)-n*t)    1      0
                        0           0  cos(n*t)
    ];

H = @(t) [
                 sin(n*t)/n       (1-cos(n*t))*2/n      0  
             (1-cos(n*t))*-2/n (4*sin(n*t)-3*n*t)/n     0
                    0                    0          sin(n*t)/n
    ];

S = @(t) [
                3*n*sin(n*t)      0         0 
              -6*n*(1-cos(n*t))   0         0
                     0            0     -n*sin(n*t)
    ];

T = @(t) [
                cos(n*t)         2*sin(n*t)         0
              -2*sin(n*t)       4*cos(n*t)-3        0
                    0                 0          cos(n*t)
    ];

%initial values
x0 = [
    100 %x
    100 %y
    10  %z
    ];

v0 = [
    2   %x'
    2   %y'
    1   %z'
    ];

%solutions
x = @(t) M(t)*x0 + H(t)*v0;
v = @(t) S(t)*x0 + T(t)*v0;

solx = x0;
solv = v0;

for i = 2:tf
    solx(:,i) = x(i);
    solv(:,i) = v(i);
end

%fEul
solFeul = [x0
           v0];
%Heun
solHeun = [x0
           v0];

for k = 1:N
 solFeul(:,k+1) = solFeul(:,k) + h*CW(solFeul(:,k),n); %fEul
 solHeun(:,k+1) = solHeun(:,k) + (h/2)*CW(solHeun(:,k),n) + (h/2)*CW(solHeun(:,k) + h*CW(solHeun(:,k),n),n); %Heun
end

%graphs
figure; set(gcf,'position',[10,10,1000,700])
%top sx plot
subplot(2,2,1);
plot3(solv(1,:),solv(2,:),solv(3,:),'k-'); grid on; hold on;
plot3(solFeul(4,:),solFeul(5,:),solFeul(6,:),'bo'); grid on; hold on;
plot3(solHeun(4,:),solHeun(5,:),solHeun(6,:),'go'); grid on; hold on;
legend('velocity path','fEul simulation','Heun simulation','location','best')
xlabel("x'"); ylabel("y'"); zlabel("z'");
title('Velocity evolution in time',sprintf('Time: %i Nodes: %i Steps size: %0.5g',tf,N,h));
%top dx plot
subplot(2,2,2);
plot3(solx(1,:),solx(2,:),solx(3,:),'k-'); grid on; hold on;
plot3(solx(1,tf),solx(2,tf),solx(3,tf),'ro'); hold on;
plot3(solFeul(1,:),solFeul(2,:),solFeul(3,:),'bo'); grid on; hold on;
plot3(solHeun(1,:),solHeun(2,:),solHeun(3,:),'go'); grid on; hold on;
legend('chaser path','end','fEul simulation','Heun simulation','location','best')
xlabel('x'); ylabel('y'); zlabel('z');
title('Position evolution in time',sprintf('Time: %i Nodes: %i Steps size: %0.5g',tf,N,h));
%bottom sx plot
subplot(2,2,3);
plot(1:tf,solv(1,:),'k-',1:tf,solv(2,:),'r-',1:tf,solv(3,:),'b-'); hold on;
plot(0:h:tf,solFeul(4,:),'y-',0:h:tf,solFeul(5,:),'g-',0:h:tf,solFeul(6,:),'c-'); hold on;
legend("x'","y'","z'","x' fEul","y' fEul","z' fEul",'location','best')
xlabel('time'); ylabel('velocity variables');
title('Evolution of velocity variables')
%bottom dx plot
subplot(2,2,4);
plot(1:tf,solx(1,:),'k-',1:tf,solx(2,:),'r-',1:tf,solx(3,:),'b-'); hold on;
plot(0:h:tf,solFeul(1,:),'y-',0:h:tf,solFeul(2,:),'g-',0:h:tf,solFeul(3,:),'c-'); hold on;
legend("x","y","z","x fEul","y fEul","z fEul",'location','best')
xlabel('time'); ylabel('position variables');
title('Evolution of position variables')


%CV equations
function dydt = CW(y,n)
        dydt = zeros(6,1);
        dydt(1) = y(4);
        dydt(2) = y(5);
        dydt(3) = y(6);
        dydt(4) = 3*n^2*y(1) + 2*n*y(5);
        dydt(5) = -2*n*y(4);
        dydt(6) = -1*n^2*y(3);
end



