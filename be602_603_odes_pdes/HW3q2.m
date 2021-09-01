clear all; 
close all;
clc;

%Boundary Stuff
a = 1; %x bound
b = 2; %y bound 
xdomain = -a:0.05:a;
ydomain = 0:0.05:b;
[Xmesh, Ymesh] = meshgrid(xdomain, ydomain);
Lindex = 0:1:5;
Mindex = 0:1:5;

%Array initializations
C = zeros(length(Lindex),length(Mindex)); %preallocated for speed
syms x;
PLexps = @(x) ... %define array of callable Legendre polys 
    [0.*x + 1;                                          %l = 0
    x;                                                  %l = 1
    .5.*(3.*x.^2 - 1);                                  %l = 2
    .5.*(5.*x.^3 - 3.*x);                               %l = 3
    .125.*(35.*x.^4 - 30.*x.^2 + 3);                    %l = 4
    .125.*(63.*x.^5 - 70.*x.^3 + 15.*x)];               %l = 5

Pindex = tril(ones(length(Lindex),length(Mindex)));
PL = Pindex .* PLexps(x); %gives lower tri matrix of Pl polys

PLMesh = PLexps(Xmesh); %gives crazy matrix of concantanated mesh values

%FS coeff and FS summation
f = 0;
for Lcount = 1:1:length(Lindex)
    for Mcount = 1:1:length(Mindex)
        L = Lindex(Lcount);
        M = Mindex(Mcount);
        if M == 0
            C(Lcount,Mcount) = ((2*L+1)/4) * int(PL(Lcount,Mcount),0,a);
        else
            C(Lcount,Mcount) = ((2*L+1)/(M*pi))*sin(M*pi/2)*int(PL(Lcount,Mcount),0,a);
        end
        f = f + ((C(Lcount,Mcount) .* PLMesh((length(Xmesh)*L+1):(length(Xmesh)*L+length(Xmesh)),:) .* cos((M*pi/b).*Ymesh)));
    end
end

C(abs(C) < 1E-14) = 0; %removes weird values that are essentially zero

for i = 1:length(C(:,1))
    for j = 1:length(C(1,:))
        disp({'Fourier Coefficients: ', string(C(i,j))}) %echos the coeffs to the command window
    end
end

%Surf Plot
figure();
fxyhandle = surf(Xmesh, Ymesh, f); %this surf plot does not obey the Dirichlet boundaries, but idk what my error is
colormap('jet');
xlabel('x');
ylabel('y');
title('Surf plot of f (Legendre)');
caxis([-.5 1.5])
surf_colorbarhandle = colorbar;
set(get(surf_colorbarhandle, 'YLabel'),'String', 'Value of f');

%Pcolor Plot
figure();
pcolorhandle = pcolor(Xmesh, Ymesh, f); 
colormap('jet');
xlabel('x');
ylabel('y');
title('pcolor plot of f (Legendre)');
caxis([-.5 1.5])
pcolor_colorbarhandle = colorbar;
set(get(pcolor_colorbarhandle, 'YLabel'),'String', 'Value of f');









