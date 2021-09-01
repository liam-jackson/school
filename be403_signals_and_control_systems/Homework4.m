%Liam Jackson
%HW 4

%2e

n = 0:20;
x = ones(21,1);
h = ones(21,1);

x([11:21],1) = 0;
h([11:21],1) = 0;
H = zeros(length(n));
H(:,1) = h;

for i = 2:length(n)
    Hi = H(:,i-1);
    Hi = circshift(Hi,1);
    Hi(1) = 0;
    H(:,i) = Hi;
end

y = .1*H*x;
circshift(y,1);
y(1) = 0;
stem(n,y)
ylim([-.5,1.5])
grid on
title('Discretized Convolution of x(t) and h(t)')

%3

clear all;

figure
subplot(4,1,1)
load ('handel.mat');
x = y;
%sound(x)
plot(x)
title('Handel')

subplot(4,1,2)
load('capture1.mat');
%sound(h)
plot(h)
title('Book in a Box (Response)')

subplot(4,1,3)
y = conv(x,h);
%sound(y)
plot(y)
title('Handel in a Box')

clear x
subplot(4,1,4)
load('laughter.mat');
x = y;
z = conv(x,h);
%sound(z)
plot(z)
title('Laughter in a Box');
%this terrifying clip shows what multiple people in a large plywood box
%would sound like while laughing maniacally

















