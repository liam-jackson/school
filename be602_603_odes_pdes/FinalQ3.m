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

