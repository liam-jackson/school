close all, clear all, clc;
warning('off','MATLAB:polyfit:RepeatedPointsOrRescale')

%% Part 1

%%% Importing / Sorting Data
[x1, x2, y] = textread('besseldata.txt', ' %f%f%f', 'headerlines', 1);
data_raw = sortrows([sqrt(x1.^2 + x2.^2), y]);

r = data_raw(:, 1);
y = data_raw(:, 2);

%%% Bessel Approx
k_bes = 1;
bes_approx = besselj(0, k_bes*r);

fig1 = figure(1);
dot_sz = 0.2;
line_w = 2.5;
scatter(r, y, dot_sz, '.');
title({'Timpanic Memb Displacement', 'approximated by Bessel Fxn (J_0)'});
hold on;
plot(r, bes_approx, 'LineWidth', line_w);
hold off;
xlabel('r');
ylabel('Intensity');
legend({'Real Data', 'J_0'});
%%
%%% Build Matrix of Polyfit coefficients, polyval values 
max_poly_order = 14;
n_coeffs = max_poly_order + 1;

coeffs_arr = zeros([n_coeffs, max_poly_order]); %15x14
poly_vals = zeros([length(r), max_poly_order]); %5000x14

for poly_order_ind = 1:max_poly_order
    poly_coeffs = polyfit(r, y, poly_order_ind);    
    for r_ind = 1:length(poly_coeffs)
        coeffs_arr(r_ind, poly_order_ind) = poly_coeffs(r_ind); %Array is in DESCENDING ORDER of poly coeffs 
    end
    poly_vals(:, poly_order_ind) = polyval(poly_coeffs, r);    
end

data_poly_vals = [r, y, poly_vals];

data_var_names = {'r', 'y'};
poly_var_nums = 1:length(poly_vals(1,:));
poly_var_names = "p = " + poly_var_nums;
var_names = {[data_var_names, poly_var_names]};

data_table = array2table(data_poly_vals,...
    'VariableNames',var_names{1});


%%
%%% Calculate Residuals
residuals_arr = zeros([max_poly_order, 2]);
for poly_order_ind = 1:max_poly_order
    y_real = data_poly_vals(:, 2);
    y_model = data_poly_vals(:, poly_order_ind + 2);
    residuals_arr(poly_order_ind, 1) = poly_order_ind;
    residuals_arr(poly_order_ind, 2) = residuals(y_real, y_model);
end

res_table = array2table(residuals_arr,...
    'VariableNames',{'Polynomial Order', 'Residual Sum'})
%%
%%% Plotting LS Poly fits
data_labels = data_table.Properties.VariableNames;

fig2 = figure(2);
sgtitle({'Membrane Displacement Data', 'vs. OLS Polynomial Fits'});
number_of_plots = max_poly_order;
r = data_table.('r');
y = data_table.('y');

for plot_id = 1:number_of_plots
    subplot(number_of_plots / 2, 2, plot_id);
    scatter(r, y, dot_sz, '.');
    hold on;
    plot(r, data_table.(string(data_labels(plot_id + 2))), 'LineWidth', line_w)
    hold off;
    legend({'y real', string(data_labels(plot_id + 2))}) 
end

%% Part 2

%%% (k = 5)-fold CV 
k_cv = 5; 
max_round_number = 20;

cumulative_res_vals = zeros([max_poly_order, max_round_number + 1]); %14x21
cumulative_res_vals(:,1) = 1:1:max_poly_order;
cumulative_coeffs = zeros([n_coeffs, max_poly_order, max_round_number]); %15x14x20

%%% 20 rnds of CV
for cv_round = 1:max_round_number
    [cumulative_res_vals(:, cv_round + 1), cumulative_coeffs(:, :, cv_round)] = res_calc(data_raw, max_poly_order, k_cv);
end

PE = [cumulative_res_vals(:, 1), cumulative_res_vals(:, 2:end) ./ (length(data_raw(:,1)) / k_cv)];

PE_table = array2table(PE,...
    'VariableNames',...
    {'Poly Ord', 'Round 1', 'Round 2', 'Round 3', 'Round 4', 'Round 5',...
    'Round 6', 'Round 7', 'Round 8', 'Round 9', 'Round 10',...
    'Round 11', 'Round 12', 'Round 13', 'Round 14', 'Round 15',...
    'Round 16', 'Round 17', 'Round 18', 'Round 19', 'Round 20'})

fig3 = figure(3);
for cv_round = 1:max_round_number
    plot(PE(:, 1), PE(:, cv_round + 1));
    hold on;
end
title('Predictive Error vs. Polynomial Model Order');
xlabel('Order');
ylabel('PE');
grid on;
hold off;

%% Part 3

% A polynomial LS-fit of order 10 seems to have the best compromise of
% accuracy and economy of variables. A substantial reduction in error
% occurs from order-9 to order-10, with no substantial decrease with
% additional (11, 12, 13, 14) order terms. 

%% Part 4

opt_ord = 10;
deg_10_coeffs = mean(squeeze(cumulative_coeffs(1:opt_ord+1, opt_ord, :)), 2) 
opt_model = polyval(deg_10_coeffs, r);

x1 = sort(x1);
x2 = sort(x2);
[x1_new, x2_new] = meshgrid(-15:.2:15);
r_new = sqrt(x1_new.^2 + x2_new.^2);
z = polyval(deg_10_coeffs, r_new);

[x1_raw, x2_raw, y_raw] = textread('besseldata.txt', ' %f%f%f', 'headerlines', 1);

model_surf = surf(x1_new, x2_new, z, 'EdgeColor', 'none');
hold on;
scatter3(x1_raw, x2_raw, y_raw, 15, '.')

%%

function [residuals, coeffs_arr] = res_calc(data_raw, max_poly_order, k_cv)
r = data_raw(:, 1);
y = data_raw(:, 2);

train_data_indices = sort(randsample(length(r), (k_cv - 1) * length(r) / k_cv)); 
removed_data_indices = true(1, length(r));
removed_data_indices(train_data_indices) = false;

train_data = data_raw(train_data_indices, :);
test_data = data_raw(removed_data_indices, :);

n_coeffs = max_poly_order + 1;

coeffs_arr = zeros([n_coeffs, max_poly_order]); %15x14
poly_vals = zeros([length(train_data(:,1)), max_poly_order]); %5000x14

for poly_order_ind = 1:max_poly_order
    poly_coeffs = polyfit(train_data(:, 1), train_data(:, 2), poly_order_ind);    
    for r_ind = 1:length(poly_coeffs)
        coeffs_arr(r_ind, poly_order_ind) = poly_coeffs(r_ind); %Array is in DESCENDING ORDER of poly coeffs 
    end
    poly_vals(:, poly_order_ind) = polyval(poly_coeffs, train_data(:, 1));    
end

data_poly_vals = [train_data(:, 1), train_data(:, 2), poly_vals];

%%% Calculate Residuals
residuals = zeros([max_poly_order, 2]);
for poly_order_ind = 1:max_poly_order
    y_real = data_poly_vals(:, 2);
    y_model = data_poly_vals(:, poly_order_ind + 2);
    res = abs(y_real - y_model).^2;
    residuals(poly_order_ind, 1) = poly_order_ind;
    residuals(poly_order_ind, 2) = sum(res);
end

residuals = residuals(:,2);

end

function binned_data = bin_this_data(data, k)

num_data_pts = size(data, 1);
rows_per_bin = floor(num_data_pts / k);
extra_rows_needed = mod(num_data_pts, k);

perm_ind = randperm(num_data_pts);
perm_data = data(perm_ind, :);

rowDist = rows_per_bin * ones(1, k);

for extra_row = 1:extra_rows_needed
    rowDist(extra_row) = rowDist(extra_row) + 1;
end

binned_data = mat2cell(perm_data, rowDist);

end














































