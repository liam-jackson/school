close all; clear all; clc;

x1 = linspace(0,15,16);
x2 = linspace(0,15,16);

[X1, X2] = meshgrid(x1,x2);

H = zeros(size(X1));
bnd1 = 4.5;
bnd2 = 5.5;
bnd3 = 2.5;
a1 = .5 * log(4);
a2 = .5 * log(7);
a3 = .5 * log(22/6);

figure();
hold on;
for i = 1:length(x2)
    for j = 1:length(x1)
        if x1(j) > bnd1
            h1 = 1; 
        else
            h1 = -1;
        end
        if x2(i) > bnd2
            h2 = 1; 
        else
            h2 = -1;
        end
        if x2(i) > bnd3
            h3 = 1; 
        else
            h3 = -1;
        end
        H(i,j) = sign(a1*h1 + a2*h2 + a3*h3);
        if H(i,j) == -1
            scatter(x1(j),x2(i), 100, 'rx')
        elseif H(i,j) == 1
            scatter(x1(j),x2(i), 100, 'b.')
        end
    end
end
hold off;
title('Adaboost H(x) for weak classifiers h_1, h_2, h_3');
xlabel('x1');
ylabel('x2');

figure();
surf(H);




