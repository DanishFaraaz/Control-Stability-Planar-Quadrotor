clear all
close all
clc

g = 9.81;
m = 1;
Ixx = 8.1e-3;

A = [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; 0 0 -g 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
B = [0 0; 0 0; 0 0; 0 0; 1/m 0; 0 1/Ixx];
C = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
D = 0;

sys = ss(A,B,C,D);

Q = C'*C;
R = 1;
K = lqr(A,B,Q,R);

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
title('Step Response with LQR Control')