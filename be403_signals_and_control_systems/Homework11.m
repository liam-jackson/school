%Liam Jackson
%Homework 11

t =  0:0.1:5;

%ic
Hi = tf([1],[1 1 0]);
figure(1)
impulse(Hi,t)
legend()

%iic
Hii = tf([5],[1 6 5 0]);
figure(2)
impulse(Hii,t)
legend()

%iiic
Hiii = tf([1],[1 0 -1]);
figure(3)
impulse(Hiii,t)
legend()

%ivc
Hiv = tf([1 1 5],[1 1 4 4]);
figure(4)
impulse(Hiv,t)
legend()
















