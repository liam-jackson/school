close all;
clear all;
clc;

tspan = [0 60];

M0 = 1;
Mp0 = 0;
AM0 = Mp0;
AMp0 = Mp0;

IC = [AM0, AMp0, M0, Mp0];

[t1, ODEsoln1] = ode45(@dAMpMfun1, tspan, IC);

pfrac1 = [ODEsoln1(:,2)+ODEsoln1(:,4)];
stress1 = [ODEsoln1(:,1)+ODEsoln1(:,2)];

[t2, ODEsoln2] = ode45(@dAMpMfun2, tspan, IC);

pfrac2 = [ODEsoln2(:,2)+ODEsoln2(:,4)];
stress2 = [ODEsoln2(:,1)+ODEsoln2(:,2)];

figure()
plot(t1, pfrac1);
hold on
plot(t2, pfrac2);
title({'1a. Effect of Static vs. Dynamic Reaction Constants',...
    'on Fraction of Phosphorylated Components'});
legend('Constant K1 = K6 = 1', 'Dynamic K1 = K6 = 1,0.2','Location', 'southeast')
xlabel('Time (s)');
ylabel('Fraction of Phosphorylated Components');
hold off

figure()
plot(t1,stress1); 
hold on 
plot(t2,stress2);
title({'1a. Effect of Static vs. Dynamic Reaction Constants',...
    'on Stress'});
legend('Stress with constant K1 = K6 = 1', 'Stress with dynamic K1 = K6 = 1,0.2','Location', 'southeast')
xlabel('Time (s)');
ylabel('Stress');
hold off

%%%%%

[t3, ODEsoln3] = ode45(@dAMpMfun3, tspan, IC);

pfrac3 = [ODEsoln3(:,2)+ODEsoln3(:,4)];
stress3 = [ODEsoln3(:,1)+ODEsoln3(:,2)];

[t4, ODEsoln4] = ode45(@dAMpMfun4, tspan, IC);

pfrac4 = [ODEsoln4(:,2)+ODEsoln4(:,4)];
stress4 = [ODEsoln4(:,1)+ODEsoln4(:,2)];

figure()
plot(t3, pfrac3);
hold on
plot(t4, pfrac4);
title({'1b. Effect of Static vs. Dynamic Reaction Constants',...
    'on Fraction of Phosphorylated Components',...
    '(with Additional K8 Coefficient)'});
legend('Constant K1 = K6 = 1', 'Dynamic K1 = K6 = 1,0.2','Location', 'southeast')
xlabel('Time (s)');
ylabel('Fraction of Phosphorylated Components');
hold off

figure()
plot(t3,stress3); 
hold on 
plot(t4,stress4);
title({'1b. Effect of Static vs. Dynamic Reaction Constants',...
    'on Stress',...
    '(with Additional K8 Coefficient)'});
legend('Stress with constant K1 = K6 = 1', 'Stress with dynamic K1 = K6 = 1,0.2','Location', 'southeast')
xlabel('Time (s)');
ylabel('Stress');
hold off

function dydt = dAMpMfun1(t, y)
    K1 = 1;
    K2 = .5;
    K3 = .4;
    K4 = (K3 / 4);
    K5 = K2;
    K6 = K1;
    K7 = .01;
    K8 = 0;
    
    %AM = y1, AMp = y2, M = y3, Mp = y4
    dydt = [...
    (K5*y(2)) + (K8*y(3)) - ((K6+K7)*y(1));     %AM
    (K3*y(4)) + (K6*y(1)) - ((K4+K5)*y(2));     %AMp
    (K2*y(4)) + (K7*y(1)) - ((K1+K8)*y(3))      %M
    (K4*y(2)) + (K1*y(3)) - ((K2+K3)*y(4))];    %Mp  
end

function dydt = dAMpMfun2(t, y)
    K2 = .5;
    K3 = .4;
    K4 = (K3 / 4);
    K5 = K2;
    K7 = .01;
    K8 = 0;
    
    if t <= 5
        K1 = 1;
        K6 = K1;
        %AM = y1, AMp = y2, M = y3, Mp = y4
        dydt = [...
        (K5*y(2)) + (K8*y(3)) - ((K6+K7)*y(1));     %AM
        (K3*y(4)) + (K6*y(1)) - ((K4+K5)*y(2));     %AMp
        (K2*y(4)) + (K7*y(1)) - ((K1+K8)*y(3))      %M
        (K4*y(2)) + (K1*y(3)) - ((K2+K3)*y(4))];    %Mp  
    else
        K1 = .2; 
        K6 = K1;
        %AM = y1, AMp = y2, M = y3, Mp = y4
        dydt = [...
        (K5*y(2)) + (K8*y(3)) - ((K6+K7)*y(1));     %AM
        (K3*y(4)) + (K6*y(1)) - ((K4+K5)*y(2));     %AMp
        (K2*y(4)) + (K7*y(1)) - ((K1+K8)*y(3))      %M
        (K4*y(2)) + (K1*y(3)) - ((K2+K3)*y(4))];    %Mp  
    end
end

function dydt = dAMpMfun3(t, y)
    K1 = 1;
    K2 = .5;
    K3 = .4;
    K4 = (K3 / 4);
    K5 = K2;
    K6 = K1;
    K7 = .01;
    K8 = .4;
    
    %AM = y1, AMp = y2, M = y3, Mp = y4
    dydt = [...
    (K5*y(2)) + (K8*y(3)) - ((K6+K7)*y(1));     %AM
    (K3*y(4)) + (K6*y(1)) - ((K4+K5)*y(2));     %AMp
    (K2*y(4)) + (K7*y(1)) - ((K1+K8)*y(3))      %M
    (K4*y(2)) + (K1*y(3)) - ((K2+K3)*y(4))];    %Mp        
end

function dydt = dAMpMfun4(t, y)
    K2 = .5;
    K3 = .4;
    K4 = (K3 / 4);
    K5 = K2;
    K7 = .01;
    K8 = .4;
    
    if t <= 5
        K1 = 1;
        K6 = K1;
            %AM = y1, AMp = y2, M = y3, Mp = y4
    dydt = [...
    (K5*y(2)) + (K8*y(3)) - ((K6+K7)*y(1));     %AM
    (K3*y(4)) + (K6*y(1)) - ((K4+K5)*y(2));     %AMp
    (K2*y(4)) + (K7*y(1)) - ((K1+K8)*y(3))      %M
    (K4*y(2)) + (K1*y(3)) - ((K2+K3)*y(4))];    %Mp  
    else
        K1 = .2; 
        K6 = K1;            
    dydt = [...
    (K5*y(2)) + (K8*y(3)) - ((K6+K7)*y(1));     %AM
    (K3*y(4)) + (K6*y(1)) - ((K4+K5)*y(2));     %AMp
    (K2*y(4)) + (K7*y(1)) - ((K1+K8)*y(3))      %M
    (K4*y(2)) + (K1*y(3)) - ((K2+K3)*y(4))];    %Mp  
    end
end








