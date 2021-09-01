%Liam Jackson, HW5q1 numeric version
close all;
clear all;
clc;

%Parameters for Quantitative Data
maxFSterm = 60; %defines the number of FS terms calculated
h0 = 3e-3;      %M (normalized from mM)
h1 = 0;         %M
h2 = 5e-3;      %M
a = 100e-9;     %m (normalized from nM)
b0 = 50e-9;     %m
eps = 10e-9;    %m
D = 5e-10;      %m^2/sec (normalized from cm^2/sec)

Cmatrix = zeros(1,60);
for nindex = 1:60
    n = nindex;
    big = ((h1-h2)/(n*pi)^2);
    sml = (a*h0/n/pi);
    Cmatrix(nindex) = (2/a)...
        *(big*(a*sin(n*pi*(b0-eps)/a)-n*pi*(b0-eps)*cos(n*pi*(b0-eps)/a))...
        +big*(a*sin(n*pi*(b0+eps)/a)-n*pi*(b0+eps)*cos(n*pi*(b0+eps)/a)...
            -a*sin(n*pi*(b0-eps)/a)+n*pi*(b0-eps)*cos(n*pi*(b0-eps)/a))...
        -sml*(cos(n*pi*(b0+eps)/a)-cos(n*pi*(b0-eps)/a))...
        +big*(a*sin(n*pi*a/a)-n*pi*a*cos(n*pi*a/a)...
            -a*sin(n*pi*(b0+eps)/a)+n*pi*(b0+eps)*cos(n*pi*(b0+eps)/a)));
end

%Spatial Domain
xstep = .5e-9;  %m
xspan = 0:xstep:a;

%Time Domain
tstep = 5e-9;   %seconds (normalized nanoseconds)
tfinal = 1e-6;
tspan = 0:tstep:tfinal;

%Meshgrid of Space-Time Vectors
[xgrid, tgrid] = meshgrid(xspan,tspan);

FSsum = 0;
for nindex = 1:maxFSterm
    n = nindex;
    Cn = Cmatrix(nindex);
    FSsum = FSsum + (Cn.*exp(-D.*(n*pi/a).^2.*tgrid).*sin(n.*pi.*xgrid./a));
end

r = FSsum;
s = xgrid*(h2-h1)/a;
u = r+s;

%Part 2, FS reconstruction vs defined IC
%Define h(x)
h = zeros(1,length(xspan));
for hindex = 1:length(xspan)
    xforh = xspan(hindex);
    if xforh > b0-eps && xforh < b0+eps
        h(hindex) = h0;
    else
        h(hindex) = 0;
    end
end

%Plot
fig1 = figure(1);
    plot(xspan, h)
    hold on;
    FSsumIC = u(1,:);
    plot(xspan,FSsumIC);
    hold off;
    grid on;
    title('FS Reconstruction of IC @ t=0');
    legend('IC','FS');
    xlabel('Spatial (m)');
    ylabel('Concentration (M)');
    xlim([0 a])
    ylim([-1e-3 2*h2])

%Part 3, Surf Plot for Diffusion Profile
fig2 = figure(2);
    surf(xgrid,tgrid,u,'EdgeColor','none')
    grid on;
    title('Steady State Profile')
    xlabel('Spatial (m)')
    ylabel('Time (sec)')
    zlabel('Concentration (M)')
    view([180 0])

%The SS profile hasn't reached the linear shape yet because not enough time
%has passed. The time coeffs in the exponential term demonstrate the
%dependency on time. Its effect on the solution decreases as time
%continues, but still has a significant impact by t=1usec. Increasing the
%timespan shows the linearization, which happens by approximately 5usec.
%
%The effective time constant is ((a/n/pi)^2/D). If we consider that low n
%values have the greatest impact on the curve shape, then for n=1, the time
%constant expression evaluates to 2.02e-6 seconds, which is where the
%exponential term decreases to 36.8% of its original value. If we sum all
%of the values of tau for n = 1:60, the time becomes 3.29e-6 seconds.
%Rerunning the simulation with this time span, we can show that the SS
%profile looks much more linear at 3.29e-6 seconds. Increasing the tspan
%past this point only enhances the linearity. 

%Overhead "pcolor" Plot 
fig3 = figure(3);
    surf(xgrid,tgrid,u,'EdgeColor','none')
    grid on;
    title('Concentration vs Space vs Time (pcolor-ish)')
    xlabel('Spatial (m)')
    ylabel('Time (sec)')
    zlabel('Concentration (M)')
    view([0 90])

%Overall Surf Plot
fig4 = figure(4);
    surf(xgrid,tgrid,u,'EdgeColor','none')
    grid on;
    title('Concentration Profile')
    xlabel('Spatial (m)')
    ylabel('Time (sec)')
    zlabel('Concentration (M)')

%Just Verifying
fig5 = figure(5);
    sgtitle('Comparison of u,s,r solutions')
fig5sub1 = subplot(3,1,1);
    plot(xspan, s(1,:));
    grid on;
    title('s(x)')
fig5sub2 = subplot(3,1,2);
    plot(xspan, r(1,:));
    grid on;
    title('r(x,t=0)')
    ylabel('Concentration (M)')
fig5sub3 = subplot(3,1,3);
    plot(xspan, u(1,:));
    grid on;
    title('u(x,t=0)')
    xlabel('Space (m)')
