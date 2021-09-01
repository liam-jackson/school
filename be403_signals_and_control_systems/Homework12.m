%Liam Jackson
%HW12

%1b
k1 = -5;
k2 = .51;
H1 = tf([k1],[1, .5-k1]);
H2 = tf([k2],[1, .5-k2]);

figure()
step(H1)
title('1b Stable')
legend('Stable')
figure()
step(H2)
title('1b Unstable')
legend('Unstable')

%1c
k3 = 5;
k4 = -.51;
H3 = tf([k3],[1, .5+k3]);
H4 = tf([k4],[1, .5+k4]);

figure()
step(H3)
title('1c Stable')
legend('Stable')
figure()
step(H4)
title('1c Unstable')
legend('Unstable')

%3b
clear all;
Pb = tf([2],[1 .2 1]);
figure()
step(Pb)
title('3b')

%3d
clear all
kd = 11.8/2;
kp = 4;
Pd = tf([2*kd 2*kp],[1 (.2+2*kd) (1+2*kp)]);
figure()
step(Pd)
title('3d')

%3e
kd = 11.8/2;
kp = 4;
ki = 10;

Pe = tf([0 11.8 8 20],[1 12.1 9 20]);
figure()
step(Pe)
title('3e')






