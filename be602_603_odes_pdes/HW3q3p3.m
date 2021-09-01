% \\\\\\\\\\\\\\\\\  Plotting routine \\\\\\\\\\\\\\\\\\\\\\\\\\\     
% -- Synthesize spherical objects (these will be fed into the slice 
%    function)   
[sxx, syy, szz] = sphere(360);               % -- Fine mesh sphere 
[shittyxx, shittyyy, shittyzz] = sphere(60); % -- Rough mesh sphere

%Part 1
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 
% -- Figure 1:  The surface profile of the planet ! 
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 
figure;   
my_sphere_slicehandle = slice(X2, Y2, Z2, T2, sxx, syy, szz); 
set(my_sphere_slicehandle, 'EdgeColor', 'none'); 
title('Spherical coord plot (at r = 1) of temperature at the planet surface !'); 
xlabel('x-axis'); 
ylabel('y-axis'); 
zlabel('z-axis'); 
my_colorbar_spherehandle = colorbar; 
set(get(my_colorbar_spherehandle, 'YLabel'), 'String', 'Temperature (C)'); 
caxis([-60 60]);   
% -- Limits the temperature color near the center of the planet 
%caxis([min(min(min(T)))  max(max(max(T)))])     

% -- Figure 2: The surface profile outside of the planet
% (Exterior Dirichlet region)
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% -- Plot temperature in terms of semi-tranparent slices
figure;

 %my_regular_slicehandle = slice(X, Y, Z, T, [0.5 0.8], [0 0.5], 0);
my_regular_slicehandle = slice(X2, Y2, Z2, T2, [0], [0], 0);

set(my_regular_slicehandle, 'EdgeColor', 'none', 'FaceAlpha', 0.9); 

hold on;

% -- Since we're plotting the exterior Dirichlet solution, we want to 
% mask the regions inside r = 1 with a solid (fine mesh) sphere !

my_sphere_slicehandle2 = slice(X2, Y2, Z2, T2, sxx, syy, szz);

set(my_sphere_slicehandle2, 'EdgeColor', 'none');
% -- Finally, we're gonna plot a rough mesh sphere outline over it!
my_sphere_slicehandle3 = slice(X2, Y2, Z2, T2, shittyxx, shittyyy, shittyzz); set(my_sphere_slicehandle3, 'EdgeColor', 'black', 'FaceColor', 'none');
xlabel('x-axis');
ylabel('y-axis');
zlabel('z-axis');
title('Rectangular slice plot of temperature');
my_colorbar_regularhandle = colorbar;

set(get(my_colorbar_regularhandle, 'YLabel'), 'String', 'Temperature (C)'); 
caxis([-60 60]); % -- Limits the temperature color near the center of the planet 
%caxis([0 max(max(max(T)))])

%Part 3
%Plotting theoretical h(phi) along greenwich meridan
figure;
plot(acos(xp3),T3)
hold on

for i = 1:length(xp3)
    if xp3(i) < -.5
        y(i) = 0;
    elseif xp3(i) < .5
        y(i) = 50;
    else
        y(i) = 0;
    end
end

plot(acos(xp3),y)
title('The \theta component of the planet surface temperature on \phi = 0 plane')
xlabel('\theta');
ylabel('h(\phi) and the Legendre reconstruction');
legend('theoretical h(\phi)','Legendre series approximation');
grid on
hold off 

%Part 5 
%Temperature profile outside of r = 1 sphere

% -- Plot temperature in terms of semi-tranparent slices
figure;
my_regular_slicehandle = slice(X5,Y5,Z5,T5,0,0,0); 
set(my_regular_slicehandle, 'EdgeColor', 'none', 'FaceAlpha', 0.5); 
hold on;

% -- Since we're plotting the exterior Dirichlet solution, we want to % mask the regions inside r = 1 with a solid (fine mesh) sphere !
my_sphere_slicehandle2 = slice(X5, Y5, Z5, T5, sxx, syy, szz);
set(my_sphere_slicehandle2, 'EdgeColor', 'none');

% -- plot a rough mesh sphere outline over it
% my_sphere_slicehandle3 = slice(X5, Y5, Z5, T5, shittyxx, shittyyy, shittyzz);
% set(my_sphere_slicehandle3, 'EdgeColor', 'black', 'FaceColor', 'none');
xlabel('x-axis');
ylabel('y-axis');
zlabel('z-axis');
title('Rectangular slice plot of temperature');
my_colorbar_regularhandle = colorbar;
 set(get(my_colorbar_regularhandle, 'YLabel'), 'String', 'Temperature (C)'); 
 caxis([-60 60]);

