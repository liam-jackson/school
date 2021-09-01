clear all;
close all;
clc;

%Q2
%Parameters
L = 50;         %cm (RDir BC)
tau = 70e+2;    %kg*cm/sec^2 (string tension)
% tau = 70;       %kg*m/sec^2 (string tension)
sig = .36e-5;   %kg/cm (string mass per unit length)
b0 = 10;        %cm (pulse center)
eps = 2.5;      %cm (half pulse width)
u0 = 0.1;       %cm (IC pulse height)
nmax = 60;      %number of FS coeffs to calculate

c = sqrt(tau/sig);  %Wave Propogation Speed (cm/sec)

xmax = 50;      %cm (furthest distance hand from bridge)
tfinal = 10e-3; %sec 

xstep = .1;                  %cm
xgrid = 0:xstep:xmax;       %0:1:50 cm
tstep = 0.01e-3;            %.01e-3 sec = 10e-6 = 10 usec 
tgrid = tstep:tstep:tfinal; %starts from tstep, not 0 (avoids div by 0)

[X, T] = meshgrid(xgrid,tgrid);

R = @(n) (4*u0/n/pi)*sin(n*pi*b0/L)*sin(n*pi*eps/L);
FSsum = zeros(1,nmax);
for n = 1:1:nmax
    FSsum(n) = R(n);
end

u = zeros(size(X));

for n = 1:1:nmax
    for xind = 1:1:length(xgrid)
        for tind = 1:1:length(tgrid)
            x = xgrid(xind);
            t = tgrid(tind);
            Rf = FSsum(n);
            u(tind,xind) = u(tind,xind)+(Rf*sin(n*pi*x/L)*cos(n*pi*c*t/L));
        end
    end
end

figure();
StringPos = surf(X,T,u,'EdgeColor','none');
title('String Position u(x,t)');
xlabel('Space (cm)');
ylabel('Time (sec)');
zlabel('String Position');
colormap 'jet';
% view([0 0]);        %IC profile
% view([90 0]);       %Position vs Time
% view([0 90]);       %pcolor-ish

figure();
StringIC = plot(xgrid,u(1,:));
title('Initial String Position (t=0+)');
xlabel('Space (cm)');
ylabel('String Position (cm)');

nh = 1;     %nth harmonic
f = (nh/2/L)*sqrt(tau/sig);
fprintf('The frequency of the note is: \n%.2f Hz\n',f)
%This is A440, the standard pitch for tuning an instrument. There's
%actually a theory that A432 sounds better, but I don't think so. Check out
%Adam Neely on Youtube if interested, he's got a video about it
