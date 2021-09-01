% Problem 2

clear all
% close all

%make meshgrids
X1 = 0:12;
X2 = 0:12;

[X1_mesh, X2_mesh] = meshgrid(X1,X2);

% code your weak classifying lines into the meshgrids
h1 = sign(0.5+ -(X1_mesh <=4.5));
h2 = sign(0.5+ -(X2_mesh <=5.5));
h3 = sign(0.5+ -(X2_mesh <=2.5));

alf1 = 0.5*log(4);
alf2 = 0.5*log(7);
alf3 = 0.5*log(11/3);

%Using the weights (alf# figure out what points on the plot are (+) vs (-))
H = sign(alf1*h1 + alf2*h2 + alf3*h3);

figure(1)
xlabel('x_1')
ylabel('x_2')
title('Plot of class regions determined by in-place adaboost', 'C_0 is red and C_1 is blue')
hold on
for ii = 1:length(X1)
    for jj = 1:length(X2)
        if H(ii,jj) == -1
            plot(jj-1,ii-1, 'rx', 'MarkerSize', 10, 'LineWidth', 2)
        elseif H(ii,jj) == 1
            plot(jj-1,ii-1, 'b.', 'MarkerSize', 10, 'LineWidth', 2)
        end
    end
end
hold off


% % % figure(3)
% % % my_pcolor = pcolor(H)
