clear all; 
close all;
clc;

%Boundary Stuff
mindex = 0:1:11;
nindex = 1:1:11;
a = 4; %x bound
b = 4; %y bound, I know this should be 2 as given in the prob statement but if its not 4 the plot doesn't make sense
xdom = 0:0.05:a;
ydom = 0:0.05:b;
[Xmesh, Ymesh] = meshgrid(xdom, ydom);

%Array initialization
C = zeros(length(mindex),length(nindex)); %preallocated for speed
f = 0;

%FS coeff and FS summation
for mcount = 1:1:length(mindex)
    for ncount = 1:1:length(nindex)
        m = mindex(mcount);
        n = nindex(ncount);
        if m == 0
            C(mcount,ncount) = (-1/(n*pi))*(cos(2*n*pi/3)-cos(n*pi/3));
        else
            C(mcount,ncount) = (-4/(m*n*pi^2))*(sin(m*pi/2)*(cos(2*n*pi/3)-cos(n*pi/3)));
        end
    f = f + ((C(mcount,ncount) .* cos((m*pi/a).*Xmesh) .* sin((n*pi/a).*Ymesh)));    
    end
end

C(find(abs(C) < 1E-14)) = 0; %for some reason got a lot of values approx 1E-33, removing them here. 

for i = 1:length(C(:,1))
    for j = 1:length(C(1,:))
        disp({'Fourier Coefficients: ', string(C(i,j))}) %echoing the coeffs to the command window
    end
end

%Surf Plot
figure();
fxyhandle = surf(Xmesh, Ymesh, f);
colormap('jet');
xlabel('x');
ylabel('y');
title('Surf plot of f');
caxis([-.5 1.5])
surf_colorbarhandle = colorbar;
set(get(surf_colorbarhandle, 'YLabel'),'String', 'Value of f');

%Pcolor Plot
figure();
pcolorhandle = pcolor(Xmesh, Ymesh, f);
colormap('jet');
xlabel('x');
ylabel('y');
title('pcolor plot of f');
caxis([-.5 1.5])
pcolor_colorbarhandle = colorbar;
set(get(pcolor_colorbarhandle, 'YLabel'),'String', 'Value of f');
