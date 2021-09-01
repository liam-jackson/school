%Liam Jackson, HW5 Calculations 

clear all;
clc;
%Note: I derived all expressions on paper, this file only is to document my
%calculations

% Q1
disp('Q1:');
%Group 1:
G1.Qa = 234;    %Afferent Vol Flow (nL/min)
G1.ff = .31;    %Filtration Fraction 
G1.Pg = 68;     %Glomerular Pressure (mmHg)
G1.Pb = 18;     %Bowman's Capsule Pressure (mmHg)
G1.Ca = 5.5;    %Afferent Protein Concentration (d/dL)

%Group 2:
G2.Qa = 242;    %Same units as above but for Group 2
G2.ff = .28;
G2.Pg = 53;
G2.Pb = 17;
G2.Ca = 5.5;

syms Cp;
Posm = 1.629*(Cp) + 0.2935*(Cp)^2;          %Osmotic Pressure (mmHg)
G1.piA = double(subs(Posm, Cp, G1.Ca));     %Afferent Osmotic Pressure 
G2.piA = double(subs(Posm, Cp, G2.Ca));     %Afferent Osmotic Pressure 

G1.Ce = G1.Ca/(1-G1.ff);    %Efferent Protein Concentration (g/dL)
G2.Ce = G2.Ca/(1-G2.ff);    %Efferent Protein Concentration (g/dL) 

G1.piE = double(subs(Posm, Cp, G1.Ce));     %Efferent Osmotic Pressure
G2.piE = double(subs(Posm, Cp, G2.Ce));     %Efferent Osmotic Pressure

G1.GFR = G1.Qa*G1.ff/60;    %GFR (nL/sec)
G2.GFR = G2.Qa*G2.ff/60;    
G1.Kf = G1.GFR/((G1.Pg-G1.Pb)-mean([G1.piE,G1.piA]));   %Filtration Coeff
G2.Kf = G2.GFR/((G2.Pg-G2.Pb)-mean([G2.piE,G2.piA]));   %(nL/sec/mmHg)

G1 = orderfields(G1)    %Compiled data for Group 1
G2 = orderfields(G2)    %Compiled data for Group 2

% Q2
clear all;

disp('Q2:');
Qu = 1500;      %Volumetric Flow of Urine (mL/24hrs)
Qu = Qu/1440;   %Vol Flow of Urine (mL/min)
Inu = 133;      %Conc Inulin in Urine (mg/dL)
Inu = Inu/100;  %Conc Inulin in Urine (mg/mL)
PAHu = 7;       %Conc PAH in Urine (mg/mL)
Inp = 1;        %Conc Inulin in Plasma (mg/dL)
Inp = Inp/100;  %Conc Inulin in Plasma (mg/mL)
PAHp = .015;    %Conc PAH in Plasma (mg/mL)
Crit = 0.4;     %Fraction of RBC in Blood Volume

QPAHU = Qu * PAHu               %Volume Flow of PAH Urine (mL/min)
GFR = (Inu/Inp)*Qu              %Vol Flow of Urine (mL/min)
ExPAH = Qu*PAHu                 %Excretion of PAH (mL/min)
FiltPAH = GFR*PAHp              %Filtration of PAH (mL/min)
SecPAH = ExPAH - FiltPAH        %Secretion of PAH (mL/min)

ClrPAH = (PAHu/PAHp)*Qu         %Clearance of PAH (mL/min)
RPF = ClrPAH/.90                %Renal Plasma Flow (mL/min)
RBF = RPF/(1-Crit)              %Renal Blood Flow (mL/min)
COpercentKidney = 100*RBF/4700  %Cardiac Output Fraction to Kidneys (%)

%Q3
clear all;

disp('Q3:');
Crp = 0.011e-3;                     %Conc Creatinine in Plasma (mg/L)
Cru0 = 2.3e-3;                      %Init Conc Creatinine in Urine (mg/L)
ECF = 14;                           %Extracellular Fluid Volume (L)
massMan = 80;                       %Mass of Mannitol (g)
molem = 182.17;                     %Molar Mass of Mannitol (g/mol)
Vsoln = .500;                       %Volume of Mannitol Solution (L)
Manp = (massMan/molem)/(Vsoln+ECF); %Molarity of Mannitol in Plasma (mol/L)
Cu = 1.200;                         %Osmolarity of (saturated) Urine (OsM) 
Qu0 = 0.55e-3;                      %Initial Vol Flow Rate of Urine (L/min)

GFR = (Cru0/Crp)*Qu0            %GFR (L/min)
MolExc0 = Qu0*Cu                %Initial Molar Excretion (mol/min)
ManClc = GFR*Manp               %Mannitol Clearance (mol/min)
MolExcf = MolExc0 + ManClc      %Final Molar Excretion (mol/min)
Quf = MolExcf/Cu                %Final Vol Flow Rate of Urine (L/min)
delQu = Quf-Qu0                 %Increase in Urine Flow (L/min)

%Q4
clear all;

disp('Q4:');
GFR = 120;      %mL/min
CrPro = 1.5;    %mg/min
CrSec = 0.1;    %mg/min
ECF = 14e+3;    %Volume of ECF (mL)

%a. 
Crp0 = (CrPro-CrSec)/GFR     %Plasma Concentration of Creatinine (mg/mL)

%b.
syms t;
CrPron = 1.5;       %mg/min
CrSecn = .05;       %mg/min
GFRn = 60;          %mL/min

k = (Crp0-(CrPron-CrSecn)/GFRn); %Expression from handwritten work
eq = @(t) k*exp((-GFRn/ECF)*t)+((CrPron-CrSecn)/GFRn); 
eqstring = vpa(eq,4); 
fprintf('[Cr]p = %s \n',char(eqstring)) %Print Numeric Expression 

%c. 
fprintf('[Cr]p at 30mins after blockage is: %f \n',eq(30))

