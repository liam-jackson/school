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

%%

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

%%
clear all;
close all;
clc;

%Q3a (proof of concept)
ri = 2e6;       %ohm/cm axial resistance
rm = 50e3;      %ohm/cm transmembrane resistance
cm = 150e-9;    %F/cm membrane capacitance
ro = 1e6;       %ohm/cm sheath fluid resistance

a = 0;          %cm Pulse location
b0 = 100e-3;    %sec Pulse width
Q0 = 100e-9;    %C Pulse charge

A = 1/(ri*cm);
B = 1/(rm*cm);
G = ro/((ro+ri)*cm);

% L is the leading coefficient of the solutions. 
% L = G*Q0/2/pi;                %Part a
L = G*Q0/(2*pi*sqrt(2*B));      %Part b

xstep = 0.01;       %cm dx
xmax = .25;         %cm spatial extent
tstep = 0.0005;   %sec
tfinal = 150e-3;    %sec time extent

xgrid = -xmax:xstep:xmax;
tgrid = tstep:tstep:(tfinal+b0);

[X, T] = meshgrid(xgrid,tgrid);

v = zeros(size(X'));

vgreen = @(x,t) ((L*exp(-B*t)*sqrt(2))/(sqrt(4*A*t)))*exp(-((x-a)^2)/(4*A*t));

for i = 1:1:length(xgrid)
    for j = 1:1:length(tgrid)
        x = xgrid(i);
        t = tgrid(j);
        vimpres(i,j) = vgreen(x,t)*heaviside(t);
        v(i,j) = vgreen(x,t)*heaviside(t) + vgreen(x,t-b0)*heaviside(t-b0);
    end
end

excess = b0/tstep;
X = X(1:end-excess,:);
T = T(1:end-excess,:);
v = v'; 
v = v(1:end-excess,:); 
vimpres = vimpres'; 
vimpres = vimpres(1:end-excess,:); 

figure();
ImpRes = surf(X,T,vimpres,'edgecolor','none');
title('Impulse Response');
xlabel('Spatial (cm)');
ylabel('Time (sec)');
zlabel('Voltage (V)');

%3b (I think this is what my answer should look like after convolving
%impulse response with new k_ext signal)
kext = zeros(length(X(1,:)),length(T(:,1)));
for kind = 1:1:length(T(:,1))
    ttemp = tgrid(kind);
    kext(:,kind) = Q0*(heaviside(ttemp) - heaviside(ttemp-b0));
end

%this is just to illustrate what my expectations were for the results
%of the math on paper
for i = 1:length(X(1,:))
    for j = 1:length(T(:,1))
    soln(:,i) = conv(vimpres(:,i),kext(i,:)); %doesn't count, used conv()    
    end
end

figure();
Ideal = surf(soln,'edgecolor','none');
title('Voltage across membrane with step function (Ideal)');
xlabel('Spatial (arb)'); %arbitrary units because conv() screwed the vector lengths and I don't want to fix it
ylabel('Time (arb)');
zlabel('Voltage (arb)');

figure();
Idealpc = pcolor(soln);
Idealpc.EdgeColor = 'none';
title('pcolor plot (Ideal)');
xlabel('Space (arb)');
ylabel('Time (arb)');
%%
figure();
Vsurf = surf(X,T,v,'EdgeColor','none');
title('V(x,t) with 2 impulses (part a)');
xlabel('Spatial (cm)');
ylabel('Time (sec)');
zlabel('Voltage (V)');

figure();
sgtitle('Summed impulses from part a');
fig2sub1 = subplot(2,1,1);
VIC = surf(X,T,v,'EdgeColor','none');
title('V(x,t=0)');
xlabel('Spatial (cm)');
ylabel('Time (sec)');
zlabel('Voltage (V)');
view([0 0])

fig2sub2 = subplot(2,1,2);
V_of_t = surf(X,T,v,'EdgeColor','none');
title('V(t)');
xlabel('Spatial (cm)');
ylabel('Time (sec)');
zlabel('Voltage (V)');
view([90 0])

%Below is what really happened. Obviously either my code sucks, or my math
%is broken. In any case, I feel dumb given how much time I spent trying to
%get this right

%Temp Conv Mesh
bgrid = xgrid;
% v = zeros(length(tgrid),length(xgrid));
for xind = 1:1:length(X(1,:))
    for tind = 1:1:length(T(:,1))            
        xtemp = X(tind,xind);
        ttemp = T(tind,xind);
        btemp = bgrid;
        fxn1 = exp(-sqrt(B/A)*abs(btemp));
        fxn2 = exp(-((xtemp-btemp-a).^2)./(4.*A.*ttemp));
        f1f2 = fxn1.*fxn2.*heaviside(ttemp);
        fxn3 = exp(-sqrt(B/A).*abs(btemp));
        fxn4 = exp(-((xtemp-btemp-a).^2)./(4.*A.*(ttemp-b0)));
        f3f4 = fxn3.*fxn4.*heaviside(ttemp-b0);
        combo1 = (1/sqrt(ttemp)).*trapz(bgrid,f1f2);
        combo2 = (1/sqrt(ttemp-b0)).*trapz(bgrid,f3f4);
        combot = combo1 - combo2;
        v(tind,xind) = L*exp(-B*ttemp).*combot;
    end
end

figure();
myVsurf = surf(X,T,v,'EdgeColor','none');
title('My V(x,t) Surf');
xlabel('Spatial (cm)');
ylabel('Time (sec)');
zlabel('Voltage (V)');

figure();
Idealpc = pcolor(v);
Idealpc.EdgeColor = 'none';
title('pcolor plot (My Case)');
xlabel('Space (arb)');
ylabel('Time (arb)');

figure();
sgtitle('My case (not ideal)');
fig4sub1 = subplot(2,1,1);
VIC = surf(X,T,v,'EdgeColor','none');
title('V(x,t=0)');
xlabel('Spatial (cm)');
ylabel('Time (sec)');
zlabel('Voltage (V)');
view([0 0])

fig4sub2 = subplot(2,1,2);
myV_of_t = surf(X,T,v,'EdgeColor','none');
title('V(t)');
xlabel('Spatial (cm)');
ylabel('Time (sec)');
zlabel('Voltage (V)');
view([90 0])

%That's all I got, I'm eager to know how you'd go about solving/coding
%this!! Thank you!
