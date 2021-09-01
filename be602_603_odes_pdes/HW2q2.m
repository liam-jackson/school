%Liam Jackson
%HW2 question 2

close all
clear all
% clc
 
%Time vector:
tinit = 0;
dt = .1; 
tfinal = 60;
trigid = [tinit:dt:tfinal];

%Parameters:
a = .7;
b = .8;
c = 1.9;
global constants
constants = [a, b, c];

%Initial Conditions:
x0 = 2;
y0 = -1;
IC = [x0 y0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1a.
%ODEs:
[tsol, fsolRaw] = ode45(@ODEsys1, [tinit tfinal], IC);
xsolRaw = fsolRaw(:,1);
ysolRaw = fsolRaw(:,2);

%Interpolation:
xsol = interp1(tsol,fsolRaw(:,1),trigid)';
ysol = interp1(tsol,fsolRaw(:,2),trigid)';

%Plot:
figure()
plot(xsol, ysol, 'o');
title({'Problem 1a. My chosen value is (c = 1.9)','y(t) vs. x(t)'});
xlabel('x(t)')
ylabel('y(t)')
axis square;
xmina = -2.5; 
xmaxa = 2.5; 
ymina = -2.5;
ymaxa = 2.5;
axis([xmina xmaxa ymina ymaxa]);
grid on 

%1b.
figure()
plot(trigid, xsol, 'o');
title({'Problem 1b. My chosen value is (c = 1.9)', 'x(t) vs. t'});
xlabel('t')
ylabel('x(t)')
xmina = 0; 
xmaxa = trigid(length(trigid)); 
ymina = -2;
ymaxa = 2;
axis([xmina xmaxa ymina ymaxa]);
grid on 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2.
%Parameter Update:
global constants2
c2 = 3;
constants2 = [a, b, c2];

%ODEs:
[tsol2, fsolRaw2] = ode45(@ODEsys2, [tinit tfinal], IC);
xsolRaw2 = fsolRaw2(:,1);
ysolRaw2 = fsolRaw2(:,2);

%Interpolation:
xsol2 = interp1(tsol2,fsolRaw2(:,1),trigid)';
ysol2 = interp1(tsol2,fsolRaw2(:,2),trigid)';

%Plot:
%2a.
figure()
plot(xsol2, ysol2, 'o');
title('2a. y(t) vs. x(t)');
xlabel('x(t)')
ylabel('y(t)')
axis square;
xmina = -2.5; 
xmaxa = 2.5; 
ymina = -2.5;
ymaxa = 2.5;
axis([xmina xmaxa ymina ymaxa]);
grid on 

%2b.
figure()
plot(trigid, xsol2, 'o');
title('2b. x(t) vs. t');
xlabel('t')
ylabel('x(t)')
xmina = 0; 
xmaxa = trigid(length(trigid)); 
ymina = -2;
ymaxa = 2;
axis([xmina xmaxa ymina ymaxa]);
grid on 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3.
%ODEs:
[tsol3, fsolRaw3] = ode45(@ODEsys3, [tinit tfinal], IC);
xsolRaw3 = fsolRaw3(:,1);
ysolRaw3 = fsolRaw3(:,2);

%Interpolation:
xsol3 = interp1(tsol3,fsolRaw3(:,1),trigid)';
ysol3 = interp1(tsol3,fsolRaw3(:,2),trigid)';

%Plot:
%3a.
figure()
plot(xsol3, ysol3, 'o');
title('3a. y(t) vs. x(t)');
xlabel('x(t)')
ylabel('y(t)')
axis square;
xmina = -2.5; 
xmaxa = 2.5; 
ymina = -2.5;
ymaxa = 2.5;
axis([xmina xmaxa ymina ymaxa]);
grid on 

%3b.
figure()
plot(trigid, xsol3, 'o');
title('3b. x(t) vs. t');
xlabel('t')
ylabel('x(t)')
xmina = 0; 
xmaxa = trigid(length(trigid)); 
ymina = -2;
ymaxa = 2;
axis([xmina xmaxa ymina ymaxa]);
grid on 

function fdt = ODEsys1(t, f)
    fdt = zeros(2,1);
    z1 = -.4;
    global constants;
        a = constants(1);
        b = constants(2);
        c = constants(3);

    fdt =   [c * (f(2) + f(1) - (f(1)^3)/3 + z1);
            -(1/c) * (f(1) - a + b * f(2))];        
end

function fdt2 = ODEsys2(t, f2)
    fdt2 = zeros(2,1);
    z2 = -.4;
    global constants2; 
        a2 = constants2(1);
        b2 = constants2(2);
        c2 = constants2(3);
    
    if (t >= 18 && t < 40)
        fdt2 =  [c2 * (f2(2) + f2(1) - (f2(1)^3)/3);
                -(1/c2) * (f2(1) - a2 + b2 * f2(2))];    
    else
        fdt2 =  [c2 * (f2(2) + f2(1) - (f2(1)^3)/3 + z2);
                -(1/c2) * (f2(1) - a2 + b2 * f2(2))]; 
    end
end

function fdt3 = ODEsys3(t, f3)
    fdt3 = zeros(2,1);
    z3 = -.4;
    global constants2; 
        a3 = constants2(1);
        b3 = constants2(2);
        c3 = constants2(3);
    
    if (t >= 22 && t < 40)
        fdt3 =  [c3 * (f3(2) + f3(1) - (f3(1)^3)/3);
                -(1/c3) * (f3(1) - a3 + b3 * f3(2))];    
    else
        fdt3 =  [c3 * (f3(2) + f3(1) - (f3(1)^3)/3 + z3);
                -(1/c3) * (f3(1) - a3 + b3 * f3(2))]; 
    end
end
