%Liam Jackson HW5

clear all; clc;

%Question 3d.

t = 0:.01:3;
x = sawtooth(2*pi*t);

fxn = @(n) sin(n*2*pi*t)*(2*(sin(n*pi)-(n*pi))/(n^2*pi^2));
harm1 = fxn(1);
harm2 = fxn(2);
harm3 = fxn(3);
harm4 = fxn(4);
harm5 = fxn(5);
approx = 0;
for k = 1:100
    approx = approx + fxn(k);
end

figure(1)
sgtitle({'Harmonics of a Sawtooth Wave','Fourier Series Approximation (n = [1 .. 5])'})
h1 = subplot(5,1,1);
plot(t,harm1);
ylim([-1 1])
title('Harmonic 1')
h2 = subplot(5,1,2);
plot(t,harm2);
ylim([-1 1])
title('Harmonic 2')
h3 = subplot(5,1,3);
plot(t,harm3);
ylim([-1 1])
title('Harmonic 3')
h4 = subplot(5,1,4);
plot(t,harm4);
ylim([-1 1])
title('Harmonic 4')
h5 = subplot(5,1,5);
plot(t,harm5);
ylim([-1 1])
title('Harmonic 5')
xlabel('Time(s)','FontSize',18)
p1 = get(h1,'position');
p2 = get(h2,'position');
p3 = get(h3,'position');
p4 = get(h4,'position');
p5 = get(h5,'position');
height = .5*(p1(2)-p5(2));
ydim = axes('position',[p1(1) height+.25*p5(2) 5*p1(3) p1(4)],'visible','off');
ylabel('x(t)','visible','on','FontSize',18);

figure(2)
plot(t,x)
hold on
plot(t,approx)
title({'Sawtooth Wave vs. Fourier Series Approximation','(n = [1 .. 100])'})
xlabel('Time(s)','FontSize',18)
ylabel('x(t)','visible','on','FontSize',18);
legend('Sawtooth Wave','Fourier Series')

%4
clear all;
p = load('pulsespec.mat');

freq = p.pulse_spectrum(:,1);
mag = p.pulse_spectrum(:,2);
t = 0:.001:2.5;

figure(3)
sgtitle('Fourier Series of Pulse Spectrum')
subplot(2,2,1)
synth = 0;
for i = 1
    for j = 1
        synth = synth + mag(i,1)*sin(freq(j,1)*t*2*pi);
    end
end
plot(t,synth)
title('n = [1]')
xlabel('Time (s)')
ylabel('x(t)')
grid on

subplot(2,2,2)
synth = 0;
for i = 1:3
    for j = 1:3
        synth = synth + mag(i,1)*sin(freq(j,1)*t*2*pi);
    end
end
plot(t,synth)
title('n = [1..3]')
xlabel('Time (s)')
ylabel('x(t)')
grid on

subplot(2,2,3)
synth = 0;
for i = 1:10
    for j = 1:10
        synth = synth + mag(i,1)*sin(freq(j,1)*t*2*pi);
    end
end
plot(t,synth)
title('n = [1..10]')
xlabel('Time (s)')
ylabel('x(t)')
grid on

subplot(2,2,4)
synth = 0;
for i = 1:length(mag)
    for j = 1:length(freq)
        synth = synth + mag(i,1)*sin(freq(j,1)*t*2*pi);
    end
end
plot(t,synth)
title('n = [1..351]')
xlabel('Time (s)')
ylabel('x(t)')
grid on








