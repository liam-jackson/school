%Homework 9 Liam Jackson

%4b

t = 0:.01:10;

R = .8;
C = 1.2;

x1 = @(t) cos(6*t);
x2 = @(t) sin(6*t);
h = @(t) (1/C)*exp(-t/(R*C))

yss1 = h(t).*x1(t);
yss2 = h(t).*x2(t);

plot(t,yss1)
hold on 
plot(t,yss2)
title('Cosine vs Sine Decay to Steady State')
xlabel('Time (s)')
ylabel('Voltage (V)')
hold off


