%Liam Jackson, HW5q1 symbolic version
close all;
clear all;
clc;

%Parameters for Integration
syms h0 h1 h2 a n b0 eps x t;

%Define 3 Intervals of Integration
x1dom = [0 b0-eps];
x2dom = [b0-eps b0+eps];
x3dom = [b0+eps a];

%Define Transient Solutions over each region
r1 = ((h1-h2)/a)*x;
r2 = ((h1-h2)/a)*x + h0;
r3 = ((h1-h2)/a)*x;

%Calculate Integrals for FS Coefficients
cnp1 = int(r1*sin(n*pi*x/a),x1dom);
cnp2 = int(r2*sin(n*pi*x/a),x2dom);
cnp3 = int(r3*sin(n*pi*x/a),x3dom);

%Compile Symbolic Expression for FS Coeffs
cn = (2/a)*simplify(cnp1 + cnp2 + cnp3); %included a/2 sine IP value here

%Parameters for Quantitative Data
maxFSterm = 60; %defines the number of FS terms calculated
h0 = 3e-3;      %M (normalized from mM)
h1 = 0;         %M
h2 = 5e-3;      %M
a = 100e-9;     %m (normalized from nM)
b0 = 50e-9;     %m
eps = 10e-9;    %m
D = 5e-10;      %m^2/sec (normalized from cm^2/sec)

%Spatial Domain
xstep = .5e-9;  %m
xspan = 0:xstep:a;

%Time Domain
tstep = 5e-9;   %seconds (normalized nanoseconds)
tfinal = 1e-6;
tspan = 0:tstep:tfinal;

%Meshgrid of Space-Time Vectors
[xgrid, tgrid] = meshgrid(xspan,tspan);

%Calculating FS Coefficients and 
%Fourier Series Summation for r(x,0)
Cmatrix = zeros(1,maxFSterm);
FSsum = 0;
for nindex = 1:maxFSterm
    n = nindex;
    Cmatrix(nindex) = double(subs(cn));
    Cn = Cmatrix(nindex);
    FSsum = FSsum + (Cn.*exp(-D.*(n*pi/a).^2.*tgrid).*sin(n.*pi.*xgrid./a));
end

s = (((h2-h1)/a).*xgrid + h1);
r = FSsum;
u = s + r;

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
%%
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
