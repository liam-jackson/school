clear all;
close all;
clc;

%Problem 2

b = 10;
r = 0:0.1:b;

I1 = trapz(r,r.*besselj(0,(2.40.*r./b)))
I2 = trapz(r,r.^2.*besselj(1,(3.83.*r./b)))
norm2 = (1/2).*b^2.*besselj(1,(2.40.*b./b)).^2