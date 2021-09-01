close all;
clear all;
clc;
syms R phi theta xdumb;

%Allowed value ranges for L,M 
Lallow = 0:1:6;
Mallow = 0:1:2;%max(Lallow);

%Limits of integration for functions h(phi), g(theta), g(x)
phispan = [0 2*pi];
thetaspan = [pi/3 2*pi/3];
xspan = [cos(thetaspan(2)) cos(thetaspan(1))];

%Functions after Sep of Vars of System
g = 50; %corellated to thetaspan above, which 
        %defines the range of theta (and eventually x) where g = 50
h = (1/2).*(1+cos(2.*phi));

%Initializing Legendre Poly Matrix
Pval = xdumb^0.*ones(length(Lallow),length(Mallow));

%Calculating Coefficient Matrix
c = zeros(length(Lallow), length(Mallow));
for Lindex = 1:1:max(Lallow)+1
    for Mindex = 1:1:max(Mallow)+1
        L = Lallow(Lindex);
        M = Mallow(Mindex);
        Pval(Lindex,Mindex) = ((1/((2^L)*factorial(L)))...
            *(1-xdumb^2)^(M/2)*diff((xdumb^2-1)^L,L+M));
        gPint = int(g.*Pval, xspan);
        if M <= L
            if M == 0
                %hint is the h(phi) integral
                hint = int(h,phispan); 
                c(Lindex,Mindex) = ...
                    ((2*L+1)/4/pi)*(factorial(L-M)/factorial(L+M))...
                    *hint...
                    *gPint(Lindex,Mindex);
            else
                %hint is the h(phi)*cos(M*phi) integral
                hint = int(h.*cos(M*phi),phispan); 
                c(Lindex,Mindex) = ...
                    ((2*L+1)/2/pi)*(factorial(L-M)/factorial(L+M))...
                    *hint...
                    *gPint(Lindex,Mindex);
            end
        end     
    end
end

%Coeff Matrix (Container) for Numeric Values
c = double(c); 
