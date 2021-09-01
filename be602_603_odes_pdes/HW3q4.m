clear all;
close all;
clc;

%Variables and Params
E0 = .5;
a = 1;
x = -3:.01:3;
y = 0;
z = -3:.01:3;
[X, Z] = meshgrid(x,z);
Y = y;

r = sqrt((X.^2)+(Y.^2)+(Z.^2));
Vtotal = ((-E0/2).*r.^2).*(1-((a^5).*(r.^-5))).*((3.*Z.^2)-1);

voltageContours = [-1.5:.1:-.1, -.05, 0, .05, .1:.1:1.5];

figure()
[ContourHandle, h] = contour(X, Z, Vtotal, voltageContours);
clabel(ContourHandle,h);
CfinalColorbarHandle = colorbar;
title('Equipotential Contour Plot')
xlabel('X')
ylabel('Z')
xlim([min(x) max(x)]);
zlim([min(z) max(z)]);

clear x z X Z;
E0 = .5;
a = 1;
x = -3:.2:3;
z = -3:.2:3;
[X, Z] = meshgrid(x,z);

r = sqrt(X.^2+Y.^2+Z.^2);
Vtotal = ((-E0/2).*r.^2).*(1-((a^5).*(r.^-5))).*((3.*Z.^2)-1);
delX = diff(X(1,1:2));
delZ = diff(Z(1:2,1));
[VgradX, VgradZ] = gradient(Vtotal, delX, delZ);

figure()
[ContourHandle, h] = contour(X, Z, Vtotal, voltageContours);
hold on 
quiver(X,Z,-VgradX,-VgradZ)
CfinalColorbarHandle = colorbar;
title('Equipotential lines with associated E-field vectors')
xlabel('X')
ylabel('Z')
xlim([min(x) max(x)]);
zlim([min(z) max(z)]);
hold off

