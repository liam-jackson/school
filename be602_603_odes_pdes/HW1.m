%Liam Jackson HW1

clear all, close all, clc

%1.3:
E = 10;
eta = 1;
w = linspace(0,1E6);
syms t;

AmpEpsilon = (-1 ./ (sqrt((w .* eta).^2 + E^2)));

figure;
Q1plot = loglog(w,AmpEpsilon);
title({'Max Amplitude of Epsilon' ,'Steady State Component'})
xlabel('w (rad/s)')
ylabel('Displacement')
xlim([0 1E6])
grid on

%3.1:

%Parameters:
n = 4; % # subunits
krpos = 1; % forward rxn rate constant
krneg = .8; % reverse rxn rate constant
s0 = 1; % substrate concentration
t = [0: .02: 1]; %time vector

%ODE Sys Matrix
A = [-n*krpos*s0, krneg;
    n*krpos*s0, -krneg];

%Initial Conditions:
IC = [1e-3; 0]; % @ t = 0

%Eigenvalues/corresponding eigenvectors
[V,D] = eig(A); %V vectors, D values
lambda1 = D(1,1);
lambda2 = D(2,2);
evect1 = V(:,1);
evect2 = V(:,2);

ICRowReduction = [evect1, evect2, IC]; %combining paper work with 
c = rref(ICRowReduction);                   %RREF of this matrix
c1 = c(1,3);
c2 = c(2,3);

%Answer:
R = c1 .* evect1 * exp(lambda1 * t) + c2 .* evect2 * exp(lambda2 *t);

R0 = R(1,:);
R1 = R(2,:);

%3.2:
figure; 
R0vsT = plot(t, R0);
hold on
R1vsT = plot(t, R1);
legend('R0','R1','Location', 'best');
xmin = 0;
xmax = 1;
ymin = 0;
ymax = 1E-3;
axis([xmin xmax ymin ymax]);
title({'Concentrations of','R0 (Free Enzyme) and R1 (Complex)','vs time'})
xlabel('time')
ylabel('Concentration of Enzyme Form')
grid on
hold off

%3.3:
figure; 
R1vsR0 = plot(R0,R1,'o');
hold on 

%3.4:
line([0, evect1(1,1)],[0, evect1(2,1)],'Linestyle','--', 'Color', 'magenta')
line([0, evect2(1,1)],[0, evect2(2,1)],'Linestyle','--','Color', 'green')

%3.5:

%B
ICB = [.4E-3;0];
ICRowReductionB = [evect1, evect2, ICB]; %combining paper work with 
cB = rref(ICRowReductionB);                   %RREF of this matrix
c1B = cB(1,3);
c2B = cB(2,3);

RB = c1B .* evect1 * exp(lambda1 * t) + c2B .* evect2 * exp(lambda2 *t);

R0B = RB(1,:);
R1B = RB(2,:);
R1BvsR0B = plot(R0B,R1B,'+');

%C
ICC = [.2E-3;0];
ICRowReductionC = [evect1, evect2, ICC]; %combining paper work with 
cC = rref(ICRowReductionC);                   %RREF of this matrix
c1C = cC(1,3);
c2C = cC(2,3);

RC = c1C .* evect1 * exp(lambda1 * t) + c2C .* evect2 * exp(lambda2 *t);

R0C = RC(1,:);
R1C = RC(2,:);
R1CvsR0C = plot(R0C,R1C,'.');

%D
ICD = [0;.6E-3];
ICRowReductionD = [evect1, evect2, ICD]; %combining paper work with 
cD = rref(ICRowReductionD);                   %RREF of this matrix
c1D = cD(1,3);
c2D = cD(2,3);

RD = c1D .* evect1 * exp(lambda1 * t) + c2D .* evect2 * exp(lambda2 *t);

R0D = RD(1,:);
R1D = RD(2,:);
R1DvsR0D = plot(R0D,R1D,'square');

title({'Eigenvectors dictate','asymptotic behavior','of ODE system'})
legend('IC','Eigenvector 1','Eigenvector 2','ICB','ICC','ICD','Location','southeast');
axis square;
xlabel('R0');
ylabel('R1');
grid on
xmin = -1E-3;
xmax = 1E-3;
ymin = -1E-3;
ymax = 1E-3;
axis([xmin xmax ymin ymax]);

%4.1  - I really had trouble here, and resorted to recycling your code
%example from blackboard. It didn't seem legal but its the only way I got
%anything to work

clear all;

global A;   

% -- Define A-matrix

%Parameters
m1 = 1;
m2 = 1; 
k1 = 1; 
k2 = 2; 
k3 = 1; 
b1 = 0; 
b2 = 1;
b3 = 0;
f1 = 0;
f2 = 0;

%ODE system
x1dot = [];
x2dot = [];
y1dot = [];
y2dot = [];
wdot = [x1dot;x2dot;y1dot;y2dot];

x1 = [];
x2 = [];
y1 = [];
y2 = [];
w = [x1;x2;y1;y2];

g = [0;0;f1/m1;f2/m2];

A = [0,0,1,0;
    0,0,0,1;
    -(k1+k2)/m1, k2/m1, -(b1+b2)/m1, b2/m1;
    k2/m2, -(k2+k3)/m2, b2/m2, -(b2+b3)/m2];

% -- Define times

t_initial = 0;
t_final = 10*2*pi;
deltaT = 0.01*2*pi;   % -- We will utilize this in the "interp1" command
                %    so that we'll have equally-spaced time pts !
                
% -- ICs
x10 = 1;
y10 = 0;
x20 = 0;
y20 = 0;
w0 = [x10;y10;x20;y20];

% -- Solve ODE:    u'  =   Au
[t, w] = ode45('hello2', [t_initial  t_final],  w0);

% -- Set up equally-spaced time mesh for the phase portrait !

t_grid = t_initial : deltaT : t_final;  % -- Equal time points

w1_grid = interp1(t, w(:,1), t_grid);   % -- Interpolate u1(t) for equal time points
w2_grid = interp1(t, w(:,2), t_grid);   % -- Interpolate u2(t) for equal time points  
w3_grid = interp1(t, w(:,3), t_grid);   
w4_grid = interp1(t, w(:,4), t_grid);   

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% --   Figure 1:  Plot  the phase plot
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
xmin = -1;
xmax = 1;
ymin = -1;
ymax = 1;

figure;

x2vsx1 = plot(w1_grid,  w2_grid,  '.');
set(x2vsx1, 'Color', 'red', 'Marker', 'o', 'MarkerSize', 12, 'Linewidth', 2);
grid on;

xlabel('x1(t)');
ylabel('x2(t)');
title('Phase Plot for x2(t) vs x1(t)');
axis([xmin xmax ymin ymax]);
legend()

figure;

x1vst = plot(t_grid, w1_grid);
xlabel('t');
ylabel('x1(t)');
title('x1(t) vs Time');

%I ran out of time to finish the rest of the problem...















