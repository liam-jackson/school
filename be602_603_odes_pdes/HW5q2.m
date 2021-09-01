%Liam Jackson, HW5q2
clear all; 
close all; 
clc;

%Parameters 
%   [Note: velocities stored in vector to iterate through for
%   solving various profiles below]
a = 100e-9;             %m
b0 = 40e-9;             %m
eps = 10e-9;            %m
h0 = 3e-3;              %M
h1 = 0;                 %M
h2 = 5e-3;              %M
D = 5e-10;              %m^2/sec
vcm = [0, 0.2, 1, 6];   %cm/sec
v = vcm*1e-2;           %m/sec
tfinal = 1e-6;          %seconds

xstep = 1e-9;           %m
tstep = .5e-9;          %seconds
xspan = 0:xstep:a;      
tspan = 0:tstep:tfinal;

%Define Meshgrid 
[xgrid, tgrid] = meshgrid(xspan, tspan);

%Have to initialize this here
U = zeros(length(xspan),length(tspan),length(v));

%This loop iterates through the velocities and outputs associated solutions 
for vind = 1:length(v)
    v0 = v(vind);

    %Initialize Conductances
    G0 = D/(xstep^2);
    Gflow = v0/xstep;
    Gt = 1/tstep;

    %Define IC
    hx = zeros(1,length(xspan),length(v));
    for xind = 1:1:length(xspan)
        xtemp = xspan(xind);
        if xtemp >= b0-eps && xtemp <= b0+eps
            hx(1,xind,vind) = h0;
        else
            hx(1,xind,vind) = 0;
        end
    end

    %Build G Matrix
    diag = -2*G0 + Gflow + Gt;
    diagR = G0 - Gflow;
    diagL = G0;

    Gftcs = eye(size(xgrid));
    for i = 1:length(Gftcs(:,1)+1)
        for j = 1:length(Gftcs(1,:)+1)
            if Gftcs(i,j) == 1 && j == 1
                Gftcs(i,j) = diag;
                Gftcs(i,j+1) = diagR;        
            elseif Gftcs(i,j) == 1 && j == length(Gftcs(1,:)) 
                Gftcs(i,j) = diag;
                Gftcs(i,j-1) = diagL;
            elseif Gftcs(i,j) == 1 && i > 1 && j > 1
                Gftcs(i,j) = diag;
                Gftcs(i,j-1) = diagL;
                Gftcs(i,j+1) = diagR;
            end
        end
    end

    Gftcs = eye(length(xspan)-2);
    for i = 1:length(Gftcs(:,1)+1)
        for j = 1:length(Gftcs(1,:)+1)
            if Gftcs(i,j) == 1 && j == 1
                Gftcs(i,j) = diag;
                Gftcs(i,j+1) = diagR;        
            elseif Gftcs(i,j) == 1 && j == length(Gftcs(1,:)) 
                Gftcs(i,j) = diag;
                Gftcs(i,j-1) = diagL;
            elseif Gftcs(i,j) == 1 && i > 1 && j > 1
                Gftcs(i,j) = diag;
                Gftcs(i,j-1) = diagL;
                Gftcs(i,j+1) = diagR;
            end
        end
    end

    U(:,1,vind) = hx(1,:,vind);
    U(1,2:end,vind) = h1;
    U(end,2:end,vind) = h2;

    bsources = zeros(length(xspan)-2,1);
    bsources(1,1) = G0*h1;
    bsources(end,1) = (G0-Gflow)*h2;

    for ii = 2:1:length(xspan)-1
        for jj = 2:1:length(tspan)
            Uold = U(2:end-1,jj-1,vind);
            U(2:end-1,jj,vind) = (Gftcs*Uold + bsources)*tstep;
        end
    end

    s = D*(tstep/(xstep^2));
    c = v0*tstep/xstep;

    if c^2 < 2*s && 2*s < 1
        fprintf('CFL Conditions satisfied for v0 = %.1f cm/sec\n',vcm(vind))
    else
        fprintf('Change your xstep or tstep, CFL invalidfor v0 = %.1f cm/sec\n',vcm(vind));
    end
end

%Separate out the individual solutions
U1 = U(:,:,1)';
U2 = U(:,:,2)';
U3 = U(:,:,3)';
U4 = U(:,:,4)';

%Plotting Routine
fig1 = figure(1);
    colormap jet;
    pc = pcolor(xgrid,tgrid,U4);
    set(pc, 'EdgeColor', 'none');
    title('Pcolor plot for 1D Diffusion with Advection (v0 = 6cm/sec)');
    xlabel('Space (m)');
    ylabel('Time (sec)');
    c = colorbar;
    ylabel(c,'Concentration (M)');

fig2 = figure(2);
sgtitle({'Pcolor Concentration Profiles for Various','Advection Velocities'});
sub1fig2 = subplot(2,2,1);
    pc = pcolor(xgrid,tgrid,U1);
    set(pc, 'EdgeColor', 'none');
    title('No Advection (v0 = 0cm/sec)');
    xlabel('Space (m)');
    ylabel('Time (sec)');
    c = colorbar;
    ylabel(c,'Concentration (M)');
sub2fig2 = subplot(2,2,2);
    pc = pcolor(xgrid,tgrid,U2);
    set(pc, 'EdgeColor', 'none');
    title('(v0 = 0.2cm/sec)');
    xlabel('Space (m)');
    ylabel('Time (sec)');
    c = colorbar;
    ylabel(c,'Concentration (M)');
sub3fig2 = subplot(2,2,3);
    pc = pcolor(xgrid,tgrid,U3);
    set(pc, 'EdgeColor', 'none');
    title('(v0 = 1cm/sec)');
    xlabel('Space (m)');
    ylabel('Time (sec)');
    c = colorbar;
    ylabel(c,'Concentration (M)');
sub4fig2 = subplot(2,2,4);
    pc = pcolor(xgrid,tgrid,U4);
    set(pc, 'EdgeColor', 'none');
    title('(v0 = 6cm/sec)');
    xlabel('Space (m)');
    ylabel('Time (sec)');
    c = colorbar;
    ylabel(c,'Concentration (M)');

fig3 = figure(3);
    SS1 = U1(end,:);
    SS2 = U2(end,:);
    SS3 = U3(end,:);
    SS4 = U4(end,:);
    ssprof = plot(xspan,SS1,xspan,SS2,xspan,SS3,xspan,SS4);
    title({'Steady State Concentration Profiles for','Various Initial Velocities'});
    legend('v = 0cm/s','v = 0.2cm/s','v = 1cm/s','v = 6cm/s','Location','northwest');
    xlabel('Space (m)');
    ylabel('Concentration (M)');
