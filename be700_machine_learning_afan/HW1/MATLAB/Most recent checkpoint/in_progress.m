close all, clear all, clc;
warning('off','MATLAB:polyfit:RepeatedPointsOrRescale')

%%% Part 1

%%% Importing / Sorting Data
[x1, x2, y] = textread('besseldata.txt', ' %f%f%f', 'headerlines', 1);

data_xxy = sortrows([x1, x2, y]); 
data_ry = sortrows([sqrt(x1.^2 + x2.^2), y]);
data_rnormy = [normalize(data_ry(:, 1)), y];

r = data_ry(:, 1);
r_norm = data_rnormy(:, 1);
y = data_ry(:, 2);

%%% Bessel Approx
k_bes = 1;
bes_approx = besselj(0, k_bes*r);

dot_sz = 0.2;
line_w = 2.5;

% fig1 = figure(1);
% scatter(r, y, dot_sz, '.');
% hold on;
% plot(r, bes_approx, 'LineWidth', line_w);
% hold off;
% title({'Timpanic Memb Displacement', 'approximated by Bessel Fxn (J_0)'});
% xlabel('r');
% ylabel('Intensity');
% legend({'Real Data', 'J_0'});

%%% Polynomial Approximations
max_poly_order = 14;

data_poly_struct = iter_poly_vals(r, y, max_poly_order);
data_poly_arr = data_poly_struct.data_poly_vals_arr;
data_poly_table = data_poly_struct.data_poly_vals_table;

%%% Calculate Residuals
res_table = res_table(data_poly_arr)

%%% Plotting LS Poly fits
data_labels = data_poly_table.Properties.VariableNames;

fig2 = figure(2);
sgtitle({'Membrane Displacement Data', 'vs. OLS Polynomial Fits'});
number_of_plots = max_poly_order;
r = data_poly_table.('r');
y = data_poly_table.('y');

for plot_id = 1:number_of_plots
    subplot(number_of_plots / 2, 2, plot_id);
    scatter(r, y, dot_sz, '.');
    hold on;
    plot(r, data_poly_table.(string(data_labels(plot_id + 2))), 'LineWidth', line_w)
    hold off;
    legend({'y real', string(data_labels(plot_id + 2))}) 
end
%%
%%% Part 2

cv_rounds = 20;
k_cv = 5; 
for round = 1:cv_rounds
    binned_data = bin_this_data(data_ry, k_cv);
    bin_indices = 1:1:k_cv;
    for cv_round = 1:k_cv
        test_data = cell2mat(binned_data(cv_round));
        train_bins = bin_indices(1:end ~= cv_round);
        train_data = cell2mat(binned_data(train_bins));
        
%         figure(cv_round);
%         scatter(train_data(:,1), train_data(:,2)); 
%         hold on; 
%         scatter(test_data(:,1), test_data(:,2)); 
%         hold off;
        
%          data_poly_vals_cv 

    end
end














