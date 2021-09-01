close all;
clear all;
clc;

%1.
%Parameters
ri = 2e+6;      %Ohms/cm
ro = 0;         %Ohms/cm
rm = 10e+3;     %Ohms/cm
cm = 150e-9;    %Farads/cm

a = 2.5;    %cm
b = 1;      %cm
D = 1/(ri*cm); %axial charge migration
B = 1/(rm*cm); %transmemb charge leakage

%Grid Variables
xstep = a/1000;     %cm
xgrid = -a:xstep:a;

tfinal = 5e-3;      %sec
tstep = tfinal/100; %sec
tgrid = tstep:tstep:tfinal;

[X, T] = meshgrid(xgrid,tgrid);

%Define h(xgrid) double pulse
h = zeros(1,length(xgrid));
for hind = 1:1:length(xgrid)
    xtemp = xgrid(hind);
    if abs(xtemp) > b
        h(hind) = 0;
    elseif xtemp >= -b && xtemp < 0
        h(hind) = -1;
    elseif xtemp >= 0 && xtemp <= b 
        h(hind) = 1;
    end
end

%Convolving h with f
%Temp Conv Mesh
bgrid = xgrid;
soln = zeros(length(tgrid),length(xgrid));
for xind = 1:1:length(xgrid)
    for tind = 1:1:length(tgrid)            
        xtemp = X(tind,xind);
        ttemp = T(tind,xind);
        btemp = bgrid;
        f = h(xind).*exp((-(xtemp-btemp).^2)./(4*D*ttemp));    
        soln(tind,xind) = (1./(sqrt(4*pi*D.*ttemp)))*exp(-B*ttemp).*trapz(bgrid,f);
        
    end
end

%Plotting:
figure();
VoltageSurf = surf(X,T,soln,'EdgeColor','none');
hold on

%This block adds some time lines to make the plot a little more defined for
%checking variances with time. Comment out if it sucks
timelinespacing = 10;
for i = 1:timelinespacing:length(tgrid) %Plots a line every 10 data points in time dimension
    TimeLines = plot3(X(i,:),T(i,:),soln(i,:),'-k');
end

xlabel('Spatial (cm)');
ylabel('Time (sec)');
zlabel('Voltage (V)');
title({'Transmembrane Voltage (Axon Unit Length)','Response to External Double Pulse Stimulus'});
colormap jet;
c = colorbar;
set(get(c,'title'),'string','Voltage (V)');
% view([0 0]); %Uncomment this line to view the V(x) profile 
% view([90 0]); %Uncomment this line to view the V(t) profile
%Note: you can only use one of those lines at a time

