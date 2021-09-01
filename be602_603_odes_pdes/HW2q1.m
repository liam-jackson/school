%Liam Jackson
%HW2, part 1

clear all

%Toroidal mesh parameters
R0 = 2;     %major radius
r0 = 1;     %minor radius
xmina = -3;  %axis shit
xmaxa = 3;
ymina = -3; 
ymaxa = 3;
zmina = -3;
zmaxa = 3;

%Azimuth/Polar anglular domain
theta1domain = 0: 2*pi/100 : 2*pi;
theta2domain = 0: 2*pi/100 : 2*pi;

%Azimuth and polar angle grid 
[theta1grid, theta2grid] = ...
    meshgrid(theta1domain,theta2domain);

%Torus grid coordinates
X=(R0+r0.*cos(theta2grid)).*cos(theta1grid);
Y=(R0+r0.*cos(theta2grid)).*sin(theta1grid);
Z=r0.*sin(theta2grid);

%1a 
% ---------- Parameters: ----------
%Time vector:
tinitial = 0;
tfinal = 20;
dt = 0.02; 
tgrid = (tinitial : dt : tfinal)';

%Constants:
K1 = 5.75; 
w1 = 2*pi;
K2 = 0;
w2 = .2*pi;
global constants
constants = [w1, K1, w2, K2];

%Initial Conditions:
theta10 = 0;    %theta1(t=0) = 0 
theta20 = pi;   %theta2(t=0) = pi 
IC1 = [theta10, theta20];

%ODE system solution
[tsol, thetaRaw] = ode45(@ODEsys,[tinitial tfinal],IC1);

%Interpolated solution
theta1sol = interp1(tsol,thetaRaw(:,1),tgrid);
theta2sol = interp1(tsol,thetaRaw(:,2),tgrid);
thetasol = [theta1sol theta2sol];

%Torus solution coordinates
Xsol = (R0 + r0 .* cos(theta2sol)) .* cos(theta1sol);
Ysol = (R0 + r0 .* cos(theta2sol)) .* sin(theta1sol);
Zsol = r0 .* sin(theta2sol);
SolCoords = [Xsol Ysol Zsol];

%IC point
XsolInit = [Xsol(1)];
YsolInit = [Ysol(1)];
ZsolInit = [Zsol(1)];
ICpt = [XsolInit;
    YsolInit;
    ZsolInit];

%Final point
XsolFinal = [Xsol(length(tgrid))];
YsolFinal = [Ysol(length(tgrid))];
ZsolFinal = [Zsol(length(tgrid))];
Finalpt = [XsolFinal;
    YsolFinal;
    ZsolFinal];

%Draw torus grid and overlay ODE solution
figure()
torus_surfhandle = surf(X,Y,Z);
set(torus_surfhandle, 'FaceColor', [0.9 0.9 0.9], ...
    'EdgeColor', [.75 .75 .75], 'FaceAlpha', 0.2, ...
    'Linewidth', 0.1);
hold on
plot3(SolCoords(:,1), SolCoords(:,2), SolCoords(:,3), ...
    '.r', 'MarkerSize', 6)
plot3(Finalpt(1), Finalpt(2), Finalpt(3), '.g', 'MarkerSize', 70);
plot3(ICpt(1), ICpt(2), ICpt(3), '.b', 'MarkerSize', 70);
xlabel('x axis');
ylabel('y axis');
zlabel('z axis');
grid on;
axis([xmina xmaxa ymina ymaxa zmina zmaxa]);
hold off

FirstTwoPoints = [XsolInit, YsolInit, ZsolInit;
                  Xsol(2), Ysol(2), Zsol(2)];
Displacement1 = pdist(FirstTwoPoints);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1b. 
xminb = 0;      %axis shit
xmaxb = 2*pi;
yminb = -1; 
ymaxb = 12;

phi = [0: .01 : 2*pi];
phidot = (w1 - w2) - K1 * sin(phi);

figure()
phasediff_handle = plot(phi, phidot);
title({'1b. Phase Difference ODE',' (d/d\theta)\phi vs \phi'}) 
xlabel('\phi')
ylabel('(d/d\theta)\phi')
axis([xminb xmaxb yminb ymaxb])
grid on

%I think this section could be used to find the fixed points
syms phi
%The following returns: phi = asin((36*pi)/115)
%                       phi = pi - asin((36*pi)/115)
phidot = (w1 - w2) - K1 * sin(phi);
phizero = solve(phidot==0,phi)  
%The following displays the values of phi, 
%one of which I think is a stable fixed point                                
disp(['Zero of Phi = ', num2str(double(phizero(1)))])           
disp(['Zero of Phi = ', num2str(double(phizero(2)))])                                

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1c. 
dynDiffRad = (thetasol(:,1)-thetasol(:,2));
dynDiffDeg = dynDiffRad .* (180/pi);

figure()
title('1c. Dynamic Phase Difference')
yyaxis left
plot(tgrid, dynDiffRad)
ylabel('Angle of Phase Shift (rad)')
yyaxis right    
ylim([dynDiffDeg(1) dynDiffDeg(length(dynDiffDeg))]) 
ylabel('Angle of Phase Shift (deg)')
xline(3, '--');
xlabel('Time (sec)')
legend('Phase Difference','t = 3sec', 'Location', 'southeast')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1d.

v1 = sin(theta1sol);
v2 = sin(theta2sol);

figure()
plot(tgrid, v1, 'or');
hold on 
plot(tgrid, v2, '.b', 'MarkerSize', 5);
xline(3, '--');
title('1d. v1(t) and v2(t)')
xlabel('Time (sec)')
ylabel('Voltage (V)')
legend('v1(t)','v2(t)','t = 3sec');
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1e.
xmine = -1;
xmaxe = 1;
ymine = -1;
ymaxe = 1;

figure()
plot(v1, v2, 'o', 'MarkerSize', 5)
title({'1e. Phase Portrait','v2(t) vs. v1(t)'})
xlabel('v1(t)')
ylabel('v2(t)')
axis([xmine xmaxe ymine ymaxe])

%Definitely phase-locked after 3sec 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Liam Jackson
%HW2, part 2

clear constants

%Toroidal mesh parameters
R0 = 2;     %major radius
r0 = 1;     %minor radius
xmina = -3;  %axis shit
xmaxa = 3;
ymina = -3; 
ymaxa = 3;
zmina = -3;
zmaxa = 3;

%Azimuth/Polar anglular domain
theta1domain = 0: 2*pi/100 : 2*pi;
theta2domain = 0: 2*pi/100 : 2*pi;

%Azimuth and polar angle grid 
[theta1grid, theta2grid] = ...
    meshgrid(theta1domain,theta2domain);

%Torus grid coordinates
X=(R0+r0.*cos(theta2grid)).*cos(theta1grid);
Y=(R0+r0.*cos(theta2grid)).*sin(theta1grid);
Z=r0.*sin(theta2grid);

%2a 
% ---------- Parameters: ----------
%Time vector:
tinitial = 0;
tfinal = 20;
dt = 0.02; 
tgrid = (tinitial : dt : tfinal)';

%Constants:
K1 = 15; 
w1 = 2*pi;
K2 = 0;
w2 = .2*pi;
global constants
constants = [w1, K1, w2, K2];

%Initial Conditions:
theta10 = 0;    %theta1(t=0) = 0 
theta20 = pi;   %theta2(t=0) = pi 
IC1 = [theta10, theta20];

%ODE system solution
[tsol2, thetaRaw2] = ode45(@ODEsys,[tinitial tfinal],IC1);

%Interpolated solution
theta1sol2 = interp1(tsol2,thetaRaw2(:,1),tgrid);
theta2sol2 = interp1(tsol2,thetaRaw2(:,2),tgrid);
thetasol2 = [theta1sol2 theta2sol2];

%Torus solution coordinates
Xsol2 = (R0 + r0 .* cos(theta2sol2)) .* cos(theta1sol2);
Ysol2 = (R0 + r0 .* cos(theta2sol2)) .* sin(theta1sol2);
Zsol2 = r0 .* sin(theta2sol2);
SolCoords2 = [Xsol2 Ysol2 Zsol2];

%IC point
XsolInit2 = [Xsol2(1)];
YsolInit2 = [Ysol2(1)];
ZsolInit2 = [Zsol2(1)];
ICpt2 = [XsolInit2;
    YsolInit2;
    ZsolInit2];

%Final point
XsolFinal2 = [Xsol2(length(tgrid))];
YsolFinal2 = [Ysol2(length(tgrid))];
ZsolFinal2 = [Zsol2(length(tgrid))];
Finalpt2 = [XsolFinal2;
    YsolFinal2;
    ZsolFinal2];

%Draw torus grid and overlay ODE solution
figure()
torus_surfhandle = surf(X,Y,Z);
set(torus_surfhandle, 'FaceColor', [0.9 0.9 0.9], ...
    'EdgeColor', [.75 .75 .75], 'FaceAlpha', 0.2, ...
    'Linewidth', 0.1);
hold on
plot3(SolCoords2(:,1), SolCoords2(:,2), SolCoords2(:,3), ...
    '.r', 'MarkerSize', 6)
plot3(Finalpt2(1), Finalpt2(2), Finalpt2(3), '.g', 'MarkerSize', 70);
plot3(ICpt2(1), ICpt2(2), ICpt2(3), '.b', 'MarkerSize', 70);
xlabel('x axis');
ylabel('y axis');
zlabel('z axis');
grid on;
axis([xmina xmaxa ymina ymaxa zmina zmaxa]);
hold off

disp({'K1 = 5.75 displacement: ', num2str(Displacement1)})
FirstTwoPoints2 = [XsolInit2, YsolInit2, ZsolInit2;
                  Xsol2(2), Ysol2(2), Zsol2(2)];
Displacement2 = pdist(FirstTwoPoints2);
disp({'K1 = 15 displacement: ', num2str(Displacement2)})

%{
The solution points are further apart, indicating that the K1 gain 
increases speed of signal acquisition. 
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2b. 
xminb = 0;      %axis shit
xmaxb = 2*pi;
yminb = -1; 
ymaxb = 12;

phi = [0: .01 : 2*pi];
phidot = (w1 - w2) - K1 * sin(phi);

figure()
phasediff_handle = plot(phi, phidot);
title({'2b. Phase Difference ODE',' (d/d\theta)\phi vs \phi'}) 
xlabel('\phi')
ylabel('(d/d\theta)\phi')
axis([xminb xmaxb yminb ymaxb])
grid on

%I think this section could be used to find the fixed points
syms phi
phidot = (w1 - w2) - K1 * sin(phi);
phizero = solve(phidot==0,phi)  
disp(['Zero of Phi = ', num2str(double(phizero(1)))])           
disp(['Zero of Phi = ', num2str(double(phizero(2)))])                                

%{
As phi*, ie the phase difference, approaches zero, the achieved frequency
of the car's radio very closely matches the signal from the radio tower.
This also means the two waves get closer to being totally overlapped. 
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2c. 
dynDiffRad = (thetasol2(:,1)-thetasol2(:,2));
dynDiffDeg = dynDiffRad .* (180/pi);

figure()
title('2c. Dynamic Phase Difference')
yyaxis left
plot(tgrid, dynDiffRad)
ylabel('Angle of Phase Shift (rad)')
yyaxis right    
ylim([dynDiffDeg(1) dynDiffDeg(length(dynDiffDeg))]) 
ylabel('Angle of Phase Shift (deg)')
xline(3, '--');
xlabel('Time (sec)')
legend('Phase Difference','t = 3sec', 'Location', 'southeast')
hold off

%{
As t->inf, the K1=15 (2c.) curve oscillates slightly. This is different 
from the K1 = 5.75 (1c.) curve, which stays much more stable around the 
asymptotic value. This means that even though the curve from 2c. minimizes 
the phase difference faster, the 1c. curve achieves a higher stability. 
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2d.

v1 = sin(theta1sol2);
v2 = sin(theta2sol2);

figure()
plot(tgrid, v1, 'or');
hold on 
plot(tgrid, v2, '.b', 'MarkerSize', 5);
xline(3, '--');
title('2d. v1(t) and v2(t)')
xlabel('Time (sec)')
ylabel('Voltage (V)')
legend('v1(t)','v2(t)','t = 3sec');
hold off

%{
The two voltages associated to each K1 gain both oscillate with the same
frequency. The phase shift is more pronounced in the curves from 1c.
There is a noticeable fluctuation in density of the points in the curves
from 2c. again indicating an oscillation of the car's radio frequency 
around the radio tower's frequency.
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2e.
xmine = -1;
xmaxe = 1;
ymine = -1;
ymaxe = 1;

figure()
plot(v1, v2, 'o', 'MarkerSize', 5)
title({'2e. Phase Portrait','v2(t) vs. v1(t)'})
xlabel('v1(t)')
ylabel('v2(t)')
axis([xmine xmaxe ymine ymaxe])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function thetadt = ODEsys(t,theta)
global constants;
thetadt = zeros(2,1);
thetadt(1) = constants(1) + constants(2)*sin(theta(2)-theta(1));
thetadt(2) = constants(3); % + constants(4)*sin(theta(1)-theta(2));
end

