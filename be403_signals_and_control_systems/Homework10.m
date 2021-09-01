%Liam Jackson 
%Homework 10

%1
H1 = tf([5],[1 16 65 50]);
H2 = tf([500 5500 5000],[1 105.1 510.5 50]);

figure(1)
sgtitle('Bode Plots')
subplot(2,1,1)
bode(H1)
title('1.a')
subplot(2,1,2)
bode(H2)
title('1.b')

%2
clear H1, clear H2;

%b
figure(2)
H = tf([100 200 100],[1 .2 100.01 10]);
bode(H)
title('Bode Plot 2.b')
%d
figure(3)
t1 = 0:0.01:1e4;
t2 = 0:.01:100;
x1 = sin(.01*t1);
x2 = sin(10*t2);

subplot(2,1,1)
lsim(H, x1, t1)
title('w = .01 rad/sec')

subplot(2,1,2)
lsim(H, x2, t2)
title('w = 10 rad/sec')

%3
clear all;
load('waves.mat')
R = 3000;
C = 1/(72000*pi);
w = 12;
tauinv = (R*C)^-1;
t = 0:1:((length(waves)-1));

H = tf([w],[w (tauinv)]);

%c
figure(4)
bode(H)

%d
figure(5)
lsim(H,waves,t)




