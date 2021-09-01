clear all;
close all;
clc;

%Problem 1

%Parameters:
b = 10;

%Note: I found a function online which calculates the zeros of Bessel
%Functions, besselzero() , used only because I didn't want to manually type
%a matrix of roots. I cross-checked against the websites provided and the
%values are accurate. 

%Establish Integer Ranges for M,P and N,Q
Mrange = 0:1:2;   
Prange = 1:1:3;   %only 3 terms needed for sufficient accuracy
Nrange = 1;         
Qrange = 1:1:3;   %only 3 terms 

%Establishing Matrix of Bessel Fxn zeros:
BesZero = zeros(length(Mrange),length(Prange));
for Mindex = 1:1:length(Mrange)
    for Pindex = 1:1:length(Prange)
        M = Mrange(Mindex);
        P = Prange(Pindex);
        BesZeroTemp = besselzero(M,P);
        %BesZero has values 
        BesZero(Mindex,Pindex) = BesZeroTemp(1,Pindex);
    end
end

k0p = BesZero(1,:)./b; %k vals for zeros of J_0 (cosine terms)
k2p = BesZero(3,:)./b; %k vals for zeros of J_2 (cosine terms)
k1q = BesZero(2,:)./b; %k vals for zeros of J_1 (sine terms)
kforsum1 = [k0p; 0, 0, 0; k2p];
kforsum2 = k1q;

%Fourier Coeffs for Group 1
C = zeros(length(Mrange),length(Prange));
for Mindex = 1:1:length(Mrange)
    for Pindex = 1:1:length(Prange)
        M = Mrange(Mindex);
        P = Prange(Pindex);
        if M == 0
            k = k0p(Pindex);
            C(Mindex,Pindex) =...
                100/(b*k*besselj(1,k*b));    
        elseif M == 1
            C(Mindex,Pindex) = 0;
        else
            k = k2p(Pindex);
            C(Mindex,Pindex) =...
                -.008*(6*b^2*besselj(4,k*b)...
                -(k*b^3*besselj(5,k*b)))...
                /(k^2*(besselj(3,k*b)^2));
        end
    end
end

%Fourier Coeffs for Group 2
D = zeros(length(Nrange),length(Qrange));
for Nindex = 1:1:length(Nrange)
    for Qindex = 1:length(Qrange)
        k = k1q(Qindex);
        D(Nindex,Qindex) =...
            -.03*(4*b*besselj(3,k*b)...
            -(k*b^2*besselj(4,k*b)))...
            /(k^2*(besselj(2,k*b))^2);
    end
end

% -- Meshgrid setup 
radius = b; 
x_step = 0.1;     
z_start = 0; 
z_step = 0.1;     
z_end = 1;   

cylinder_plot_radius = radius;  % -- Default 

% \\\\\\\\\\\\\\\\\  Define rectangular meshgrid \\\\\\\\\\\\\\\\\ 
x = -radius:x_step:radius; 
y = x; 
z = z_start:z_step:z_end;   
[X,Y,Z] = meshgrid(x, y, z);       

% \\\\\\\\\\\\\  Define ghetto spherical coordinate axes \\\\\\\\\\\\\\\\\   
% -- Define radius 
R = sqrt(X.^2 + Y.^2);  % -- Positive value   

% -- Define phi (azimuthal angle)   
PHI = atan2(Y, X);   

%Compile FS Summation
Group1Sum = 0;
for Mterm = 1:1:length(Mrange)
    for Pterm = 1:1:length(Prange)
        M = Mrange(Mterm);
        P = Prange(Pterm);
        k = kforsum1(Mterm, Pterm);
        Group1Sum = Group1Sum +...
            (C(Mterm,Pterm)...
            .*besselj(M, k.*R)...
            .*cos(M.*PHI)...
            .*exp(-k.*Z));
    end
end
Group2Sum = 0;
for Nterm = 1:1:length(Nrange)
    for Qterm = 1:1:length(Qrange)
        N = Nrange(Nterm);
        Q = Qrange(Qterm);
        k = kforsum2(Nterm,Qterm);
        Group2Sum = Group2Sum +...
            (D(Nterm,Qterm)...
            .*besselj(N, k.*R)...
            .*sin(N.*PHI)...
            .*exp(-k.*Z));
    end
end

T = Group1Sum + Group2Sum;

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 
% \\\\\\\\\  Figure 1:  Full volumetric cylinder temperature plot \\\\\\\\\\   

figure;     
% -- Plot temperature in terms of cylindrical shells   
cylinder_sectors = 314; 
[cxx, cyy, czz] = cylinder(1, cylinder_sectors);       

[czzRows, czzCols] = size(czz); % -- Find the dimension of matrix "czz" 
czz(1,:) = z_start .* ones(1, czzCols); % -- Set all values in 1st row = z_start 
czz(2,:) = (z_start + z_step) .* ones(1, czzCols);  % -- Set all values in 2nd row = z_start + z_step     
 
cxx = cylinder_plot_radius .* cxx; 
cyy = cylinder_plot_radius .* cyy;     

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 
% \\\\\\\\\\\\\\\\\\  Now we're ready to plot the cylinders !!  \\\\\\\   

% -- Set counter = 1 
count = 1;   

% --  Plot the 1st cylinder from z = [0, step].  We're gonna store all of 
%     our plot handles (ie. pointers) in a giant struct called 
%     "my_cylinder" 

my_cylinder(1).structHandle = slice(X, Y, Z, T, cxx, cyy, czz); 
set(my_cylinder(count).structHandle, 'EdgeColor', 'none' , 'FaceAlpha', 0.5); 
czz = czz + z_step; % -- Increment the cylinder segment's z-position ! 
count = count + 1;  % -- Increase counter 
hold on;   

for dummy = (z_start+z_step) : z_step : (z_end - z_step)     
    my_cylinder(count).structHandle = slice(X, Y, Z, T, cxx, cyy, czz );     
    set(my_cylinder(count).structHandle, 'EdgeColor', 'none', 'FaceAlpha', 0.5 );     
    czz = czz + z_step; % -- Increment the cylinder segment's z-position !     
    count = count + 1;  % -- Increase counter 
end   

% -- Plot the 2 "endcap" slices for the cylinder.  
my_regular_slicehandle2 = slice(X, Y, Z, T, [], [], [0 1]); 
set(my_regular_slicehandle2, 'EdgeColor', 'none'); 

set(my_regular_slicehandle2(1), 'FaceAlpha', 0.9, 'EdgeColor', 'none');  
set(my_regular_slicehandle2(2), 'FaceAlpha', 0.9, 'EdgeColor', 'none'); 
title('Full Cylindrical coord temperature plot example (Variable r, z values)'); 
xlabel('x-axis'); 
ylabel('y-axis'); 
zlabel('z-axis'); 
my_colorbar_cylinderhandle = colorbar; 
set(get(my_colorbar_cylinderhandle, 'YLabel'), 'String', 'Temperature (C)'); 
caxis([0 100]);   % -- Limits the temperature color 
colormap(jet) 
hold off; 

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 
% \\\\\\\\\  Figure 2:  Phi = pi/2 cylinder cross section plot \\\\\\\\\\   

figure;     
% -- Plot temperature in terms of cylindrical shells   
cylinder_sectors = 314; 
[cxx, cyy, czz] = cylinder(1, cylinder_sectors);       

[czzRows, czzCols] = size(czz); % -- Find the dimension of matrix "czz" 
czz(1,:) = z_start .* ones(1, czzCols); % -- Set all values in 1st row = z_start 
czz(2,:) = (z_start + z_step) .* ones(1, czzCols);  % -- Set all values in 2nd row = z_start + z_step     
 
cxx = cylinder_plot_radius .* cxx; 
cyy = cylinder_plot_radius .* cyy;     

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 
% \\\\\\\\\\\\\\\\\\  Now we're ready to plot the cylinders !!  \\\\\\\   

% -- Set counter = 1 
count = 1;   

% --  Plot the 1st cylinder from z = [0, step].  We're gonna store all of 
%     our plot handles (ie. pointers) in a giant struct called 
%     "my_cylinder" 

my_cylinder(1).structHandle = slice(X, Y, Z, T, cxx, cyy, czz); 
set(my_cylinder(count).structHandle, 'EdgeColor', 'none' , 'FaceAlpha', 0); 
czz = czz + z_step; % -- Increment the cylinder segment's z-position ! 
count = count + 1;  % -- Increase counter 
hold on;   

% --  Plot the rest of the cylinder from z = [z_start + z_step........ ....z_end]. 
%     We'll do this by incrementing the entire content of "czz" by "z_step"   

for dummy = (z_start+z_step) : z_step : (z_end-z_step)     
    my_cylinder(count).structHandle = slice(X, Y, Z, T, cxx, cyy, czz );     
    set(my_cylinder(count).structHandle, 'EdgeColor', 'none', 'FaceAlpha', 0 );     
    czz = czz + z_step; % -- Increment the cylinder segment's z-position !     
    count = count + 1;  % -- Increase counter 
end   

% -- Plot the 2 "endcap" slices for the cylinder.  Then, we'll have to hide 
%    the vertical slices from view.   
my_regular_slicehandle2 = slice(X, Y, Z, T, 0, 0, [0 1]); 
set(my_regular_slicehandle2, 'EdgeColor', 'none'); 

% -- Hide the x and y slices from the plot 
set(my_regular_slicehandle2(1), 'FaceAlpha', 0.9, 'EdgeColor', 'none');  
set(my_regular_slicehandle2(2), 'FaceAlpha', 0.0, 'EdgeColor', 'none');
set(my_regular_slicehandle2(3), 'FaceAlpha', 0.9, 'EdgeColor', 'none');
set(my_regular_slicehandle2(4), 'FaceAlpha', 0.9, 'EdgeColor', 'none');
title('Cylinder cross section temp plot at \phi = \pi / 2'); 
xlabel('x-axis'); 
ylabel('y-axis'); 
zlabel('z-axis'); 
my_colorbar_cylinderhandle = colorbar; 
set(get(my_colorbar_cylinderhandle, 'YLabel'), 'String', 'Temperature (C)'); 
caxis([0 100]) 
colormap(jet) 
hold off; 

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 
% \\\\\\\\\  Figure 3:  Phi = pi/2 cylinder cross section plot \\\\\\\\\\   

x5 = 0:.01:10;

R2 = x5;

%a.
phia = 0;
Group101a = C(1,1).*besselj(0, k0p(1).*R2).*cos(0.*phia);
Group102a = C(1,2).*besselj(0, k0p(2).*R2).*cos(0.*phia);
Group103a = C(1,3).*besselj(0, k0p(3).*R2).*cos(0.*phia);
Group121a = C(3,1).*besselj(2, k2p(1).*R2).*cos(2.*phia);
Group122a = C(3,2).*besselj(2, k2p(2).*R2).*cos(2.*phia);
Group123a = C(3,3).*besselj(2, k2p(3).*R2).*cos(2.*phia);

Group211a = D(1,1).*besselj(1, k1q(1).*R2).*sin(1.*phia);
Group212a = D(1,2).*besselj(1, k1q(2).*R2).*sin(1.*phia);
Group213a = D(1,3).*besselj(1, k1q(3).*R2).*sin(1.*phia);

%b. 
phib = -pi/2;
Group101b = C(1,1).*besselj(0, k0p(1).*R2).*cos(0.*phib);
Group102b = C(1,2).*besselj(0, k0p(2).*R2).*cos(0.*phib);
Group103b = C(1,3).*besselj(0, k0p(3).*R2).*cos(0.*phib);
Group121b = C(3,1).*besselj(2, k2p(1).*R2).*cos(2.*phib);
Group122b = C(3,2).*besselj(2, k2p(2).*R2).*cos(2.*phib);
Group123b = C(3,3).*besselj(2, k2p(3).*R2).*cos(2.*phib);

Group211b = D(1,1).*besselj(1, k1q(1).*R2).*sin(1.*phib);
Group212b = D(1,2).*besselj(1, k1q(2).*R2).*sin(1.*phib);
Group213b = D(1,3).*besselj(1, k1q(3).*R2).*sin(1.*phib);

Ta = Group101a + Group102a + Group103a +...
    Group121a + Group122a + Group123a +...
    Group211a + Group212a + Group213a;
Tb = Group101b + Group102b + Group103b +...
    Group121b + Group122b + Group123b +...
    Group211b + Group212b + Group213b;

figure;
plot(x5, Ta);
title('Bessel Approx at \phi = 0');
xlabel('r');
ylabel('Temp');

figure;
plot(x5, Tb);
title('Bessel Approx at \phi = -\pi / 2');
xlabel('r');
ylabel('Temp');










