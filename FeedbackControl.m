% Feedback based Control

clear all
close all
clc

g = 9.81;
m = 1;
Ixx = 8.1e-1;

A = [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; 0 0 -g 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
B = [0 0; 0 0; 0 0; 0 0; 1/m 0; 0 1/Ixx];
C = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
D = 0;

sys = ss(A,B,C,D);

%poles = eig(A);
%P = [-1 -2 -3 -40 -50 -60];

% Percent Overshoot
OS = 10;
% Settling Time
Ts = 6;

% Damping coeff kesi
k = -log(OS/100)/sqrt((pi^2 + (log(OS/100))^2));

% Frequency
w = 4/(k*Ts);

poly = [1 2*k*w w^2];
R = roots(poly);
P = [R(1) R(2) 10*real(R(1)) 10*real(R(1))-1 10*real(R(1))-2 10*real(R(1))-3];

%P = [-1.33+1.49j -1.33-1.49j -13.33 -14 -15 -16];
K = place(A,B,P);

Ac = (A-B*K);
Bc = B;
Cc = C;
Dc = D;


x0 = [10; 10; 0; 0; 0; 0];

sys_cl = ss(Ac,Bc,Cc,Dc);
 
t = 0:0.01:10;
r(1,:) =1*zeros(size(t));
r(2,:) =1*zeros(size(t));
[y,t,x]=lsim(sys_cl,r,t,x0);
plot(t,y(:,1),'b',t,y(:,2),'g',t,y(:,3),'r')
legend('y','z','phi');
title('Step Response with Closed Loop Feedback Control')

height = y(:,1);
v = y(:,2);
phi = y(:,3);
figure(2)
for i=1:1001
    r = rectangle('Position',[v(i) height(i) 1 1]');
    xlim([-20 20])
    ylim([-20 20])
    drawnow
    if i<1001
        delete(r)
    end
end