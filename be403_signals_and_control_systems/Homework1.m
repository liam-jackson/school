%{
Liam Jackson
U40546227
Homework 1
%}

clear all, close all

%Question 1
figure 
time = linspace(-5, 5);
v = time+1;
timereverse = -time+1;
timedelay = -time-2;
amp = timedelay*2;
plot(time, v)
hold on 
plot(time, amp)
ax = gca;
xlim([-5 5])
xlabel('Time (sec)')
ax.XAxisLocation = 'origin'; 
ylim([-5 5])
ylabel('Signal')
ax.YAxisLocation = 'origin';
title('Input Signal vs Output Signal')
legend('Input Signal', 'Output Signal','Location','southeast')
grid
hold off

clear all

%Question 2

figure
sgtitle('{\it A. exampulus} Growth Rates')

t = linspace(0,4);

subplot(2,2,1);
y_growth = 0.1 + 0.5*exp(.3*t);
plot(t,y_growth)
xlabel('Time (hrs)')
ylabel('Optical Density')
hold on
y_delayed = 0.1 + 0.5*exp(.3*(t-1));
%y_delayedgrowth = piecewise(t<1, 0.15, t>1, y_delayed);
%for some reason I get an error about using a double within fxn piecewise
%from the symbolic toolbox
plot(t,y_delayed) %I should eventually switch y_delayed with y_delayedgrowth
title('Regular vs Delayed Growth')
legend({'Regular','Delayed'},'Location','northwest')
xlabel('Time (hrs)')
ylabel('Optical Density')
grid
hold off

subplot(2,2,2);
y_double = 2*(0.1 + 0.5*exp(.3*t));
plot(t,y_double)
title('Double Cuvet Length')
legend({'\it A. exampulus'},'Location','northwest')
xlabel('Time (hrs)')
ylabel('Optical Density')
grid

subplot(2,2,3);
z_another = .07*exp(.4*t);
combined = y_growth + z_another;
plot(t,y_growth)
hold on
plot(t,z_another)
plot(t,combined)
legend({'\it A. exampulus','\it A. notherus','Combined'},'Location','northwest')
title('Combined (Additive) Growth')
xlabel('Time (hrs)')
ylabel('Optical Density')
grid
hold off

subplot(2,2,4);
z_another = .07*exp(.4*t);
multiplied = y_growth .* z_another;
plot(t,y_growth)
hold on
plot(t,z_another)
plot(t,multiplied)
legend({'\it A. exampulus','\it A. notherus','Combined'},'Location','northwest')
title('Combined (Multiplicative) Growth')
xlabel('Time (hrs)')
ylabel('Optical Density')
grid
hold off

%Question 3

H = [1 0 2 1; 3 1 0 2; 0 1 1 4; 1 6 1 0];
x = [1; 4; 2; 3];
b = H*x

H1 = H(:,1);
H2 = H(:,2);
H3 = H(:,3);
H4 = H(:,4);

lincombo = x(1,1)*H1 + x(2,1)*H2 + x(3,1)*H3 + x(4,1)*H4

%Question 4

t_1 = linspace(0,2);
t_2 = linspace(0,(2.2));
t_3 = linspace(0,1);

p_1 = 3*sin(pi*(t_1));
p_2 = 3*sin(pi*t_2-.2*pi);
p_3 = 3*sin(2*pi*t_3);

figure
subplot(2,2,1);
plot(t_1,p_1)
xlim([0 2.5])
title('Mouse 1: Respiratory Pressure')
xlabel('Time (sec)')
ylabel('Pressure')
grid

subplot(2,2,2);
plot(t_2,p_2)
hold on
scatter(0.2,0)
xlim([0 2.5])
title('Mouse 2: Respiratory Pressure')
xlabel('Time (sec)')
ylabel('Pressure')
grid
hold off

subplot(2,2,3); %this mouse is exercising. They have the highest respiration rate
plot(t_3,p_3)
xlim([0 2.5])
title('Mouse 3: Respiratory Pressure')
xlabel('Time (sec)')
ylabel('Pressure')
grid

subplot(2,2,4);

t_1upperdegree = rad2deg(2); 
t_1d = linspace(0,t_1upperdegree,360); 
t_2upperdegree = rad2deg(2.2);
t_2d = linspace(0,t_2upperdegree,360); 

p_1d = 3*sind(pi*(t_1d));
p_2d = 3*sind(pi*t_2d-(.2*180));

plot(t_1d, p_1d)
hold on 
plot(t_2d, p_2d)
legend({'Mouse 1', 'Mouse 2'}, 'Location', 'southoutside')
title('Pressure vs. Phase Angle, Mouse 1 vs Mouse 2')
xlabel('Phase Angle (theta)')
ylabel('Pressure')
ylim([-3.5 3.5])
grid
hold off 

delayP1P2 = t_2upperdegree - t_1upperdegree

clear all

%Question 5

t = linspace(-1,1);

figure
sgtitle('Comparison of exponential functions and their sum')
subplot(3,1,1);
u1 = exp(-t).*sin(2*pi*t);
[yupper1,ylower1] = envelope(u1,6,'analytic');
plot(t, u1)
hold on
plot(t, yupper1, '--')
plot(t, ylower1, '--')
xlabel('Time')
ylim([-6.5 6.5])
legend({'u1','Upper Envelope','Lower Envelope'}, 'Location', 'eastoutside')
grid
hold off

subplot(3,1,2);
u2 = exp(t).*sin(2*pi*t);
[yupper2, ylower2] = envelope(u2, 6, 'analytic');
plot(t, u2)
hold on
plot(t, yupper2, '--')
plot(t, ylower2, '--')
ylim([-6.5 6.5])
xlabel('Time')
legend({'u2','Upper Envelope','Lower Envelope'}, 'Location', 'eastoutside')
grid
hold off

subplot(3,1,3);
v = u1+u2;
[yupper3, ylower3] = envelope(v, 6, 'analytic');
plot(t, v)
hold on 
plot(t, yupper3, '--')
plot(t, ylower3, '--')
ylim([-6.5 6.5])
xlabel('Time')
legend({'u1 + u2','Upper Envelope','Lower Envelope'}, 'Location', 'eastoutside')
grid
hold off

