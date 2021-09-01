clear all;
close all;
clc;
% Problem 1b

% Parameters:
global s h q0 C R;      %global vars for passing to ode45
s = 5/6;                %sec/cardiac cycle
cycles = 5;             %plot 5 cycles in problem
h = 1/3;                %length of time valves are open
q0 = 425;               %90mL stroke volume
Paort0 = 72;            %72mmHg end diastolic Aortic Pressure
Ppulm0 = 0;             %0mmHg end systolic Pulmonary Artery Pressure
C = [1.0666, 2.2666];   %Compliance (Aorta, Pulm)
R = [.900, .075];       %Peripheral Resistance (Aorta, Pulm)

%Define Time Span
tmax = cycles*s;        %tmax req'd for 's' cycles
tstep = .01;

%Parameters for ode45
tspan = 0:tstep:tmax+tstep;   
IC = [q0, Paort0, Ppulm0];           
options = odeset('MaxStep', tstep); %forces ode45 to use dt = .01 sec

%Solution Struct
sol = ode45(@dP, tspan, IC, options);

%Time Vector 
t = sol.x(1,(1:length(sol.x)-1)); %Modified Length to fit w/ Q vector below

%Q, P vals
Q = diff(sol.y(1,:))./tstep;
Paort = sol.y(2,(1:length(Q))); 
Ppulm = sol.y(3,(1:length(Q)));

%Extraneous Info (just practice)
HR = 60/s;
BloodVolPerCycle = trapz(t,Q)/cycles;
AorticPulseP = max(max(Paort))-min(min(Paort));
SystolicBP = min(min(Paort))+(max(max(Paort))-min(min(Paort)));
DiastolicBP = max(max(Paort))-(max(max(Paort))-min(min(Paort)));
fprintf(['Heart Rate: %.2f BPM',...
    '\nStroke Volume: %.2f mL',...
    '\nAortic Pulse Pressure: %.2f mmHg',... 
    '\nCalculated Systolic BP: %.2f mmHg'...
    '\nCalculated Diastolic BP: %.2f mmHg \n'],...
    HR,BloodVolPerCycle,AorticPulseP,SystolicBP,DiastolicBP);

%Combined Data Plot
fig1 = figure(1);
    linesize = 1;
    left_color = [0 .5 0];
    right_color = [.15 .15 .15];
    set(fig1,'defaultAxesColorOrder',[left_color; right_color]);
    hold on;
    yyaxis left;
    datQfig1 = plot(t, Q,'-','LineWidth',linesize,'Color',[0 .8 0]);
    ylabel('Volume Rate (cm^3 / min)');
    yyaxis right;
    datAfig1 = plot(t, Paort,'-','LineWidth',linesize,'Color',[1 0 0]);
    datPfig1 = plot(t, Ppulm,'-','LineWidth',linesize,'Color',[1 0 1]);
    ylabel('Pressure (mmHg)');
    hold off;
    grid on;
    title('Cardiac Output Volume, Aortic Pressure, Pulmonary Pressure');
    legend('Cardiac Output','Aortic Pressure',...
        'Pulmonary Artery Pressure','Location','southoutside');
    xlabel('Time (s)');
    
function [dydt]= dP(t,y)
global C R s h q0;

Caort = C(1);
Cpulm = C(2);
Raort = R(1);
Rpulm = R(2);

%Define ODE system for Q(t), Paort(t), Ppulm(t)
if mod(t,s) <= h
    y(1) = q0*sin(pi*mod(t,s)/h);
else
    y(1) = 0;
end
    y(2) = (y(1)/Caort) - y(2)./(Raort*Caort);
    y(3) = (y(1)/Cpulm) - y(3)./(Rpulm*Cpulm);   

dydt = [    y(1);       %Q (Note: need external derivative for Q vals)
            y(2);       %Paort
            y(3)    ];  %Ppulm

end            